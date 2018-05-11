/*
 *  linux/kernel/sched.c
 *
 *  Kernel scheduler and related syscalls
 *
 *  Copyright (C) 1991, 1992  Linus Torvalds
 *
 *  1996-12-23  Modified by Dave Grothe to fix bugs in semaphores and
 *              make semaphores SMP safe
 *  1998-11-19	Implemented schedule_timeout() and related stuff
 *		by Andrea Arcangeli
 *  1998-12-28  Implemented better SMP scheduling by Ingo Molnar
 */

/*
 * 'sched.c' is the main kernel file. It contains scheduling primitives
 * (sleep_on, wakeup, schedule etc) as well as a number of simple system
 * call functions (type getpid()), which just extract a field from
 * current-task
 */

#include <linux/config.h>
#include <linux/mm.h>
#include <linux/init.h>
#include <linux/smp_lock.h>
#include <linux/nmi.h>
#include <linux/interrupt.h>
#include <linux/kernel_stat.h>
#include <linux/completion.h>
#include <linux/prefetch.h>
#include <linux/compiler.h>

#include <asm/uaccess.h>
#include <asm/mmu_context.h>

#define BITMAP_SIZE (MAX_PRIO)

typedef struct runqueue runqueue_t;

struct mlfq {
	int nr_active;
	spinlock_t *lock;
	runqueue_t *rq;
	unsigned long bitmap[BITMAP_SIZE];
	list_t queue[MAX_PRIO];
};

static struct runqueue {
	int cpu;
	spinlock_t lock;
	unsigned long nr_running, nr_switches;
	task_t *curr, *idle;
	mlfq_t *p_mlfq, mlfq_instant[1];
} runqueues[NR_CPUS] __cacheline_aligned;

// For now let's deal with only 1 cpu
#define this_rq()		(runqueues + 0)
#define task_rq(p)		(runqueues + 0)
#define cpu_rq(cpu)		(runqueues + 0)
#define cpu_curr(cpu)	(runqueues[0].curr)
#define rt_task(p)		((p)->policy != SCHED_OTHER)

#define lock_task_rq(rq,p,flags)				\
do {								\
	rq = task_rq(p);					\
	spin_lock_irqsave(&rq->lock, flags);			\
} while (0)

#define unlock_task_rq(rq,p,flags)				\
	spin_unlock_irqrestore(&rq->lock, flags)

/*
 * Multilevel feedback queue operations.
 */
static inline void dequeue_task(task_t *p, mlfq_t *p_mlfq)
{
	//printk(KERN_INFO "dequeue_task (%d)\n", p->priority);

	p_mlfq->nr_active--;
	list_del_init(&p->run_list);
	if (list_empty(p_mlfq->queue + p->priority)) {
		//printk(KERN_INFO "dequeue_task list empty (%d)\n", p->priority);
		p_mlfq->bitmap[p->priority] = 1;
	}

	//printk(KERN_INFO "~dequeue_task(%d)\n", p->priority);
}

static inline void enqueue_task(task_t *p, mlfq_t *p_mlfq)
{
	//printk(KERN_INFO "enqueue_task (%d)\n", p->priority);
	list_add_tail(&p->run_list, p_mlfq->queue + p->priority);
	p_mlfq->bitmap[p->priority] = 0;
	p_mlfq->nr_active++;
	p->p_mlfq = p_mlfq;
	//printk(KERN_INFO "~enqueue_task (%d)\n", p->priority);
}

static inline void activate_task(task_t *p, runqueue_t *rq)
{
	//printk(KERN_INFO "activate_task (%d)\n", p->priority);
	enqueue_task(p, rq->p_mlfq);
	rq->nr_running++;
	//printk(KERN_INFO "~activate_task (%d)\n", p->priority);
}

static inline void deactivate_task(task_t *p, runqueue_t *rq)
{
	//printk(KERN_INFO "deactivate_task (%d)\n", p->priority);
	rq->nr_running--;
	dequeue_task(p, p->p_mlfq);
	p->p_mlfq = NULL;
	//printk(KERN_INFO "~deactivate_task (%d)\n", p->priority);
}

void do_priority_parenting(struct task_struct *h, struct task_struct *l) {
	//printk(KERN_INFO "do_priority_parenting (%d)\n", l->priority);
	runqueue_t *rq = this_rq();

	if (h->priority < MAX_PRIO) {
		deactivate_task(l, rq);

		l->old_priority = l->priority;
		l->priority = h->priority;
		l->counter = PRIO_TO_TIMESLICE(l->priority);

		activate_task(l, rq);
	}
	//printk(KERN_INFO "~do_priority_parenting (%d)\n", l->priority);
}

void undo_priority_parenting(struct task_struct *l) {
	//printk(KERN_INFO "undo_priority_parenting (%d)\n", l->priority);
	runqueue_t *rq = this_rq();

	if (l->old_priority != -1) {
		deactivate_task(l, rq);

		l->priority = l->old_priority;
		l->old_priority = -1;
		l->counter = PRIO_TO_TIMESLICE(l->priority);

		activate_task(l, rq);
	}

	//printk(KERN_INFO "~undo_priority_parenting (%d)\n", l->priority);
}

static inline void resched_task(task_t *p)
{
	//printk(KERN_INFO "resched_task (%d)\n", p->pid);
	int need_resched;

	need_resched = p->need_resched;
	wmb();
	p->need_resched = 1;
	/*if (!need_resched) {
		smp_send_reschedule(p->processor);
	}*/

	//printk(KERN_INFO "resched_task (%d)\n", p->pid);
}

static int try_to_wake_up(task_t * p, int synchronous)
{
	//printk(KERN_INFO "try_to_wake_up\n");
	unsigned long flags;
	int success = 0;
	runqueue_t *rq;

	lock_task_rq(rq, p, flags);

	p->state = TASK_RUNNING;
	if (!p->p_mlfq) {
		if (synchronous) {
			spin_lock(&this_rq()->lock);
			activate_task(p, this_rq());
			//printk(KERN_INFO "Scuccess activ\n");
			spin_unlock(&this_rq()->lock);
		} else {
			activate_task(p, rq);
			if ((rq->curr == rq->idle) || (p->priority < rq->curr->priority)) {
				//printk(KERN_INFO "Errorneous situation try_to_wake_up\n");
				resched_task(rq->curr);
			}
		}

		success = 1;
	}

	unlock_task_rq(rq, p, flags);
	//printk(KERN_INFO "~try_to_wake_up\n");
	return success;
}

inline int wake_up_process(task_t * p)
{
	//printk(KERN_INFO "wake_up_process (%d)", p->pid);
	return try_to_wake_up(p, 0);
}

void wake_up_forked_process(task_t * p)
{
	//printk(KERN_INFO "wake_up_forked_process (%d)", p->pid);
	runqueue_t *rq = this_rq();

	spin_lock_irq(&rq->lock);
	p->state = TASK_RUNNING;
	p->priority = p->priority + 1;
	if (p->priority > MAX_PRIO - 1)
		p->priority = MAX_PRIO - 1;
	activate_task(p, rq);
	spin_unlock_irq(&rq->lock);
	//printk(KERN_INFO "~wake_up_forked_process (%d)", p->pid);
}

asmlinkage void schedule_tail(task_t *prev)
{
	spin_unlock_irq(&this_rq()->lock);
}

void handle_tick_process(task_t* p) {
	runqueue_t *rq = this_rq();
	unsigned long flags;

	spin_lock_irqsave(&rq->lock, flags);
	if ((p->priority < MAX_PRIO) && (p->policy != SCHED_FIFO) && !--p->counter) {
		p->need_resched = 1;
		dequeue_task(p, rq->p_mlfq);
		if (++p->priority >= MAX_PRIO)
			p->priority = MAX_PRIO - 1;

		p->counter = PRIO_TO_TIMESLICE(p->priority);

		//printk(KERN_INFO "Sink process %d to queue %d timeslice %d", p->pid, p->priority, p->counter);

		// Queue the task at the tail of the next queue
		enqueue_task(p, rq->p_mlfq);
	}

	//printk(KERN_INFO "handle_tick_process for %d in queue %d with timeslice %d", p->pid, p->priority, p->counter);

	spin_unlock_irqrestore(&rq->lock, flags);
}

static inline int sched_find_first_zero_bit(unsigned long bitmap[BITMAP_SIZE]) {
	//printk(KERN_INFO "sched_find_first_zero_bit\n");
	int count = BITMAP_SIZE, i =0;
	for (i = 0; i < count; ++i) {
		if (bitmap[i] == 0) {
			//printk(KERN_INFO "~sched_find_first_zero_bit (%d)\n", i);
			return i;
		}
	}

	//printk(KERN_INFO "~sched_find_first_zero_bit (%d)\n", count);
	return count;
}

extern void timer_bh(void);
extern void tqueue_bh(void);
extern void immediate_bh(void);

/*
 * scheduler variables
 */

unsigned securebits = SECUREBITS_DEFAULT; /* systemwide security settings */

extern void mem_use(void);

/*
 *	Init task must be ok at boot for the ix86 as we will check its signals
 *	via the SMP irq return path.
 */
 
struct task_struct * init_tasks[NR_CPUS] = {&init_task, };

/*
 * The tasklist_lock protects the linked list of processes.
 *
 * task->alloc_lock nests inside tasklist_lock.
 */
rwlock_t tasklist_lock __cacheline_aligned = RW_LOCK_UNLOCKED;	/* outer */

static LIST_HEAD(runqueue_head);

struct kernel_stat kstat;
extern struct task_struct *child_reaper;

#ifdef CONFIG_SMP

#define idle_task(cpu) (init_tasks[cpu_number_map(cpu)])
#define can_schedule(p,cpu) \
	((p)->cpus_runnable & (p)->cpus_allowed & (1UL << cpu))

#else

#define idle_task(cpu) (&init_task)
#define can_schedule(p,cpu) (1)

#endif

void scheduling_functions_start_here(void) { }

static void process_timeout(unsigned long __data)
{
	struct task_struct * p = (struct task_struct *) __data;

	wake_up_process(p);
}

/**
 * schedule_timeout - sleep until timeout
 * @timeout: timeout value in jiffies
 *
 * Make the current task sleep until @timeout jiffies have
 * elapsed. The routine will return immediately unless
 * the current task state has been set (see set_current_state()).
 *
 * You can set the task state as follows -
 *
 * %TASK_UNINTERRUPTIBLE - at least @timeout jiffies are guaranteed to
 * pass before the routine returns. The routine will return 0
 *
 * %TASK_INTERRUPTIBLE - the routine may return early if a signal is
 * delivered to the current task. In this case the remaining time
 * in jiffies will be returned, or 0 if the timer expired in time
 *
 * The current task state is guaranteed to be TASK_RUNNING when this 
 * routine returns.
 *
 * Specifying a @timeout value of %MAX_SCHEDULE_TIMEOUT will schedule
 * the CPU away without a bound on the timeout. In this case the return
 * value will be %MAX_SCHEDULE_TIMEOUT.
 *
 * In all cases the return value is guaranteed to be non-negative.
 */
signed long schedule_timeout(signed long timeout)
{
	struct timer_list timer;
	unsigned long expire;

	switch (timeout)
	{
	case MAX_SCHEDULE_TIMEOUT:
		/*
		 * These two special cases are useful to be comfortable
		 * in the caller. Nothing more. We could take
		 * MAX_SCHEDULE_TIMEOUT from one of the negative value
		 * but I' d like to return a valid offset (>=0) to allow
		 * the caller to do everything it want with the retval.
		 */
		schedule();
		goto out;
	default:
		/*
		 * Another bit of PARANOID. Note that the retval will be
		 * 0 since no piece of kernel is supposed to do a check
		 * for a negative retval of schedule_timeout() (since it
		 * should never happens anyway). You just have the printk()
		 * that will tell you if something is gone wrong and where.
		 */
		if (timeout < 0)
		{
			printk(KERN_ERR "schedule_timeout: wrong timeout "
			       "value %lx from %p\n", timeout,
			       __builtin_return_address(0));
			current->state = TASK_RUNNING;
			goto out;
		}
	}

	expire = timeout + jiffies;

	init_timer(&timer);
	timer.expires = expire;
	timer.data = (unsigned long) current;
	timer.function = process_timeout;

	add_timer(&timer);
	schedule();
	del_timer_sync(&timer);

	timeout = expire - jiffies;

 out:
	return timeout < 0 ? 0 : timeout;
}

/*
 *  'schedule()' is the scheduler function. It's a very simple and nice
 * scheduler: it's not perfect, but certainly works for most things.
 *
 * The goto is "interesting".
 *
 *   NOTE!!  Task 0 is the 'idle' task, which gets called when no other
 * tasks can run. It can not be killed, and it cannot sleep. The 'state'
 * information in task[0] is never used.
 */
asmlinkage void schedule(void)
{
	//printk(KERN_INFO "schedule (%d)\n", current->pid);
	struct task_struct *prev, *next;
	mlfq_t *p_mlfq;
	runqueue_t *rq;
	list_t *queue;
	int idx;
	int yieldOrWait = -1; // 0 for yield and 1 for Wait (I/O)

	BUG_ON(!current->active_mm);
need_resched_back:
	prev = current;

	if (unlikely(in_interrupt())) {
		printk("Scheduling in interrupt\n");
		BUG();
	}

	release_kernel_lock(prev, smp_processor_id());
	rq = this_rq();
	spin_lock_irq(&rq->lock);

	switch (prev->state) {
		case TASK_INTERRUPTIBLE:
			//printk(KERN_INFO "TASK_INTERRUPTIBLE");
			if (signal_pending(prev)) {
				prev->state = TASK_RUNNING;
				break;
			}
		case TASK_UNINTERRUPTIBLE:
			//printk(KERN_INFO "TASK_UNINTERRUPTIBLE");
			yieldOrWait = 1;

		default: {
			//printk(KERN_INFO "default");
			deactivate_task(prev, rq);
			break;
		}

		case TASK_RUNNING: {
			//printk(KERN_INFO "TASK_RUNNING (counter = %d)", prev->counter);
			if (prev->counter > 0)
				yieldOrWait = 0;
		}
	}

	if (unlikely(!rq->nr_running)) {
		//printk(KERN_INFO "Turn to idle task");
		next = rq->idle;
		goto switch_tasks;
	}

	if (prev != rq->idle) {
		if (yieldOrWait == 0) {
			// If counter is not yet zero and state is TASK_RUNNING it implies the task
			// gave up to cpu before expiration of time-slice..so we schedule it in RR fashion.
			dequeue_task(prev, rq->p_mlfq);

			// Restore its timeslice
			prev->counter = PRIO_TO_TIMESLICE(current->priority);
			//printk(KERN_INFO "Yield process %d to queue %d timeslice %d", current->pid, current->priority, current->counter);

			enqueue_task(prev, rq->p_mlfq);
		} else if (yieldOrWait == 1) {

			// IF the task is waiting (Interruptable or uninterruptable)..
			// for I/O or something then upgrade its queue.
			// No need to make chages to queue as it is deactivated already.
			prev->priority = prev->priority - 1;
			if (prev->priority < 0) {
				prev->priority = 0;
			}

			//printk(KERN_INFO "Task going to wait for i/o");
		}

		// Do priority parenting stuff
		if ((prev->state == TASK_INTERRUPTIBLE) || (prev->state == TASK_UNINTERRUPTIBLE)) {
			if (prev->waiting_on != NULL) {
				if (prev->priority < prev->waiting_on->holder->priority){
					do_priority_parenting(prev, prev->waiting_on->holder);
				}
			}
		} else if ((prev->state == TASK_RUNNING) && (prev->waiting_on != NULL)) {
			undo_priority_parenting(prev->waiting_on->holder);
		}

	}

	p_mlfq = rq->p_mlfq;
	idx = sched_find_first_zero_bit(p_mlfq->bitmap);
	queue = p_mlfq->queue + idx;
	next = list_entry(queue->next, task_t, run_list);

switch_tasks:
	prev->need_resched = 0;

	if (unlikely(prev == next)) {
		//printk(KERN_INFO "Same process");
		goto same_process;
	}

	kstat.context_swtch++;
	/*
	 * there are 3 processes which are affected by a context switch:
	 *
	 * prev == .... ==> (last => next)
	 *
	 * It's the 'much more previous' 'prev' that is on next's stack,
	 * but prev is set to (the just run) 'last' process by switch_to().
	 * This might sound slightly confusing but makes tons of sense.
	 */
	rq->nr_switches++;
	rq->curr = next;
	next->processor = prev->processor;

	prepare_to_switch();
	{
		struct mm_struct *mm = next->mm;
		struct mm_struct *oldmm = prev->active_mm;
		if (!mm) {
			BUG_ON(next->active_mm);
			next->active_mm = oldmm;
			atomic_inc(&oldmm->mm_count);
			enter_lazy_tlb(oldmm, next, smp_processor_id());
		} else {
			BUG_ON(next->active_mm != mm);
			switch_mm(oldmm, mm, next, smp_processor_id());
		}

		if (!prev->mm) {
			prev->active_mm = NULL;
			mmdrop(oldmm);
		}
	}

	/*
	 * This just switches the register state and the
	 * stack.
	 */
	switch_to(prev, next, prev);
	barrier();
	rq = this_rq();

same_process:
	spin_unlock_irq(&rq->lock);
	reacquire_kernel_lock(current);

	if (current->need_resched)
		goto need_resched_back;

	//printk(KERN_INFO "~schedule\n");
	return;
}

/*
 * The core wakeup function.  Non-exclusive wakeups (nr_exclusive == 0) just wake everything
 * up.  If it's an exclusive wakeup (nr_exclusive == small +ve number) then we wake all the
 * non-exclusive tasks and one exclusive task.
 *
 * There are circumstances in which we can try to wake a task which has already
 * started to run but is not in state TASK_RUNNING.  try_to_wake_up() returns zero
 * in this (rare) case, and we handle it by contonuing to scan the queue.
 */
static inline void __wake_up_common (wait_queue_head_t *q, unsigned int mode,
			 	     int nr_exclusive, const int sync)
{
	struct list_head *tmp;
	struct task_struct *p;

	CHECK_MAGIC_WQHEAD(q);
	WQ_CHECK_LIST_HEAD(&q->task_list);
	
	list_for_each(tmp,&q->task_list) {
		unsigned int state;
                wait_queue_t *curr = list_entry(tmp, wait_queue_t, task_list);

		CHECK_MAGIC(curr->__magic);
		p = curr->task;
		state = p->state;
		if (state & mode) {
			WQ_NOTE_WAKER(curr);
			if (try_to_wake_up(p, sync) && (curr->flags&WQ_FLAG_EXCLUSIVE) && !--nr_exclusive)
				break;
		}
	}
}

void __wake_up(wait_queue_head_t *q, unsigned int mode, int nr)
{
	if (q) {
		unsigned long flags;
		wq_read_lock_irqsave(&q->lock, flags);
		__wake_up_common(q, mode, nr, 0);
		wq_read_unlock_irqrestore(&q->lock, flags);
	}
}

void __wake_up_sync(wait_queue_head_t *q, unsigned int mode, int nr)
{
	if (q) {
		unsigned long flags;
		wq_read_lock_irqsave(&q->lock, flags);
		__wake_up_common(q, mode, nr, 1);
		wq_read_unlock_irqrestore(&q->lock, flags);
	}
}

void complete(struct completion *x)
{
	unsigned long flags;

	spin_lock_irqsave(&x->wait.lock, flags);
	x->done++;
	__wake_up_common(&x->wait, TASK_UNINTERRUPTIBLE | TASK_INTERRUPTIBLE, 1, 0);
	spin_unlock_irqrestore(&x->wait.lock, flags);
}

void wait_for_completion(struct completion *x)
{
	spin_lock_irq(&x->wait.lock);
	if (!x->done) {
		DECLARE_WAITQUEUE(wait, current);

		wait.flags |= WQ_FLAG_EXCLUSIVE;
		__add_wait_queue_tail(&x->wait, &wait);
		do {
			__set_current_state(TASK_UNINTERRUPTIBLE);
			spin_unlock_irq(&x->wait.lock);
			schedule();
			spin_lock_irq(&x->wait.lock);
		} while (!x->done);
		__remove_wait_queue(&x->wait, &wait);
	}
	x->done--;
	spin_unlock_irq(&x->wait.lock);
}

#define	SLEEP_ON_VAR				\
	unsigned long flags;			\
	wait_queue_t wait;			\
	init_waitqueue_entry(&wait, current);

#define	SLEEP_ON_HEAD					\
	wq_write_lock_irqsave(&q->lock,flags);		\
	__add_wait_queue(q, &wait);			\
	wq_write_unlock(&q->lock);

#define	SLEEP_ON_TAIL						\
	wq_write_lock_irq(&q->lock);				\
	__remove_wait_queue(q, &wait);				\
	wq_write_unlock_irqrestore(&q->lock,flags);

void interruptible_sleep_on(wait_queue_head_t *q)
{
	SLEEP_ON_VAR

	current->state = TASK_INTERRUPTIBLE;

	SLEEP_ON_HEAD
	schedule();
	SLEEP_ON_TAIL
}

long interruptible_sleep_on_timeout(wait_queue_head_t *q, long timeout)
{
	SLEEP_ON_VAR

	current->state = TASK_INTERRUPTIBLE;

	SLEEP_ON_HEAD
	timeout = schedule_timeout(timeout);
	SLEEP_ON_TAIL

	return timeout;
}

void sleep_on(wait_queue_head_t *q)
{
	SLEEP_ON_VAR
	
	current->state = TASK_UNINTERRUPTIBLE;

	SLEEP_ON_HEAD
	schedule();
	SLEEP_ON_TAIL
}

long sleep_on_timeout(wait_queue_head_t *q, long timeout)
{
	SLEEP_ON_VAR
	
	current->state = TASK_UNINTERRUPTIBLE;

	SLEEP_ON_HEAD
	timeout = schedule_timeout(timeout);
	SLEEP_ON_TAIL

	return timeout;
}

void scheduling_functions_end_here(void) { }

#if CONFIG_SMP
/**
 * set_cpus_allowed() - change a given task's processor affinity
 * @p: task to bind
 * @new_mask: bitmask of allowed processors
 *
 * Upon return, the task is running on a legal processor.  Note the caller
 * must have a valid reference to the task: it must not exit() prematurely.
 * This call can sleep; do not hold locks on call.
 */
void set_cpus_allowed(struct task_struct *p, unsigned long new_mask)
{
	new_mask &= cpu_online_map;
	BUG_ON(!new_mask);

	p->cpus_allowed = new_mask;

	/*
	 * If the task is on a no-longer-allowed processor, we need to move
	 * it.  If the task is not current, then set need_resched and send
	 * its processor an IPI to reschedule.
	 */
	if (!(p->cpus_runnable & p->cpus_allowed)) {
		if (p != current) {
			p->need_resched = 1;
			smp_send_reschedule(p->processor);
		}
		/*
		 * Wait until we are on a legal processor.  If the task is
		 * current, then we should be on a legal processor the next
		 * time we reschedule.  Otherwise, we need to wait for the IPI.
		 */
		while (!(p->cpus_runnable & p->cpus_allowed))
			schedule();
	}
}
#endif /* CONFIG_SMP */

#ifndef __alpha__

void set_user_nice(task_t *p, long nice)
{
	unsigned long flags;
	mlfq_t *p_mlfq;
	runqueue_t *rq;

	if (p->nice == nice)
		return;

	lock_task_rq(rq, p, flags);

	p_mlfq = p->p_mlfq;
	if (p_mlfq) {
		dequeue_task(p, p_mlfq);
	}

	p->nice = nice;
	p->priority = NICE_TO_PRIO(nice);
	if (p_mlfq) {
		enqueue_task(p, p_mlfq);
		//printk(KERN_INFO "Errorneous situation set_user_nice\n");
		resched_task(rq->curr);
	}

	unlock_task_rq(rq, p, flags);
}

/*
 * This has been replaced by sys_setpriority.  Maybe it should be
 * moved into the arch dependent tree for those ports that require
 * it for backward compatibility?
 */

asmlinkage long sys_nice(int increment)
{
	long newprio;

	/*
	 *	Setpriority might change our priority at the same moment.
	 *	We don't have to worry. Conceptually one call occurs first
	 *	and we have a single winner.
	 */
	if (increment < 0) {
		if (!capable(CAP_SYS_NICE))
			return -EPERM;
		if (increment < -40)
			increment = -40;
	}
	if (increment > 40)
		increment = 40;

	newprio = current->nice + increment;
	if (newprio < -20)
		newprio = -20;
	if (newprio > 19)
		newprio = 19;

	set_user_nice(current, newprio);
	return 0;
}

#endif

static inline struct task_struct *find_process_by_pid(pid_t pid)
{
	struct task_struct *tsk = current;

	if (pid)
		tsk = find_task_by_pid(pid);
	return tsk;
}

static int setscheduler(pid_t pid, int policy, 
			struct sched_param *param)
{
	//printk(KERN_INFO "setscheduler()\n");
	struct sched_param lp;
	struct task_struct *p;
	int retval;
	unsigned long flags;
	mlfq_t *p_mlfq;
	runqueue_t *rq;

	retval = -EINVAL;
	if (!param || pid < 0)
		goto out_nounlock;

	retval = -EFAULT;
	if (copy_from_user(&lp, param, sizeof(struct sched_param)))
		goto out_nounlock;

	/*
	 * We play safe to avoid deadlocks.
	 */
	read_lock_irq(&tasklist_lock);

	p = find_process_by_pid(pid);

	retval = -ESRCH;
	if (!p)
		goto out_unlock_tasklist;

	lock_task_rq(rq,p,flags);

	if (policy < 0)
		policy = p->policy;
	else {
		retval = -EINVAL;
		if (policy != SCHED_FIFO && policy != SCHED_RR &&
				policy != SCHED_OTHER)
			goto out_unlock;
	}

	//printk(KERN_INFO "setscheduler() 1 \n");
	/*
	 * Valid priorities for SCHED_FIFO and SCHED_RR are 1..99, valid
	 * priority for SCHED_OTHER are 0...255.
	 */
	retval = -EINVAL;
	if ((policy == SCHED_FIFO || policy == SCHED_RR) && (lp.sched_priority < 0 || lp.sched_priority > 99))
		goto out_unlock;
	if ((policy == SCHED_OTHER) && (lp.sched_priority <0 || lp.sched_priority >= MAX_PRIO))
		goto out_unlock;

	//printk(KERN_INFO "setscheduler() 2 \n");

	retval = -EPERM;
	if ((policy == SCHED_FIFO || policy == SCHED_RR) && 
	    !capable(CAP_SYS_NICE))
		goto out_unlock;
	if ((current->euid != p->euid) && (current->euid != p->uid) &&
	    !capable(CAP_SYS_NICE))
		goto out_unlock;

	p_mlfq = p->p_mlfq;
	if (p_mlfq)
		deactivate_task(p, task_rq(p));

	retval = 0;
	p->policy = policy;
	p->rt_priority = lp.sched_priority;

	p->priority = lp.sched_priority;
	if (p_mlfq)
		activate_task(p, task_rq(p));

	//printk(KERN_INFO "Success in changing priority");

	current->need_resched = 1;

out_unlock:
	unlock_task_rq(rq,p,flags);
out_unlock_tasklist:
	read_unlock_irq(&tasklist_lock);

out_nounlock:
//printk(KERN_INFO "~setscheduler()\n");
	return retval;
}

asmlinkage long sys_sched_setscheduler(pid_t pid, int policy, 
				      struct sched_param *param)
{
	return setscheduler(pid, policy, param);
}

asmlinkage long sys_sched_setparam(pid_t pid, struct sched_param *param)
{
	return setscheduler(pid, -1, param);
}

asmlinkage long sys_sched_getscheduler(pid_t pid)
{
	struct task_struct *p;
	int retval;

	retval = -EINVAL;
	if (pid < 0)
		goto out_nounlock;

	retval = -ESRCH;
	read_lock(&tasklist_lock);
	p = find_process_by_pid(pid);
	if (p)
		retval = p->policy;
	read_unlock(&tasklist_lock);

out_nounlock:
	return retval;
}

asmlinkage long sys_sched_getparam(pid_t pid, struct sched_param *param)
{
	struct task_struct *p;
	struct sched_param lp;
	int retval;

	retval = -EINVAL;
	if (!param || pid < 0)
		goto out_nounlock;

	read_lock(&tasklist_lock);
	p = find_process_by_pid(pid);
	retval = -ESRCH;
	if (!p)
		goto out_unlock;
	lp.sched_priority = p->priority;
	read_unlock(&tasklist_lock);

	/*
	 * This one might sleep, we cannot do it with a spinlock held ...
	 */
	retval = copy_to_user(param, &lp, sizeof(*param)) ? -EFAULT : 0;

out_nounlock:
	return retval;

out_unlock:
	read_unlock(&tasklist_lock);
	return retval;
}

asmlinkage long sys_sched_yield(void)
{
	set_current_state(TASK_RUNNING);

	schedule();

	return 0;
}

/**
 * yield - yield the current processor to other threads.
 *
 * this is a shortcut for kernel-space yielding - it marks the
 * thread runnable and calls sys_sched_yield().
 */
void yield(void)
{
	sys_sched_yield();
}

void __cond_resched(void)
{
	set_current_state(TASK_RUNNING);
	schedule();
}

asmlinkage long sys_sched_get_priority_max(int policy)
{
	int ret = -EINVAL;

	switch (policy) {
	case SCHED_FIFO:
	case SCHED_RR:
		ret = 99;
		break;
	case SCHED_OTHER:
		ret = 0;
		break;
	}
	return ret;
}

asmlinkage long sys_sched_get_priority_min(int policy)
{
	int ret = -EINVAL;

	switch (policy) {
	case SCHED_FIFO:
	case SCHED_RR:
		ret = 1;
		break;
	case SCHED_OTHER:
		ret = 0;
	}
	return ret;
}

asmlinkage long sys_sched_rr_get_interval(pid_t pid, struct timespec *interval)
{
	struct timespec t;
	struct task_struct *p;
	int retval = -EINVAL;

	if (pid < 0)
		goto out_nounlock;

	retval = -ESRCH;
	read_lock(&tasklist_lock);
	p = find_process_by_pid(pid);
	if (p)
		jiffies_to_timespec(p->policy & SCHED_FIFO ?
				0 : PRIO_TO_TIMESLICE(p->priority), &t);
	read_unlock(&tasklist_lock);
	if (p)
		retval = copy_to_user(interval, &t, sizeof(t)) ? -EFAULT : 0;
out_nounlock:
	return retval;
}

static void show_task(struct task_struct * p)
{
	unsigned long free = 0;
	int state;
	static const char * stat_nam[] = { "R", "S", "D", "Z", "T", "W" };

	printk("%-13.13s ", p->comm);
	state = p->state ? ffz(~p->state) + 1 : 0;
	if (((unsigned) state) < sizeof(stat_nam)/sizeof(char *))
		printk(stat_nam[state]);
	else
		printk(" ");
#if (BITS_PER_LONG == 32)
	if (p == current)
		printk(" current  ");
	else
		printk(" %08lX ", thread_saved_pc(&p->thread));
#else
	if (p == current)
		printk("   current task   ");
	else
		printk(" %016lx ", thread_saved_pc(&p->thread));
#endif
	{
		unsigned long * n = (unsigned long *) (p+1);
		while (!*n)
			n++;
		free = (unsigned long) n - (unsigned long)(p+1);
	}
	printk("%5lu %5d %6d ", free, p->pid, p->p_pptr->pid);
	if (p->p_cptr)
		printk("%5d ", p->p_cptr->pid);
	else
		printk("      ");
	if (p->p_ysptr)
		printk("%7d", p->p_ysptr->pid);
	else
		printk("       ");
	if (p->p_osptr)
		printk(" %5d", p->p_osptr->pid);
	else
		printk("      ");
	if (!p->mm)
		printk(" (L-TLB)\n");
	else
		printk(" (NOTLB)\n");

	{
		extern void show_trace_task(struct task_struct *tsk);
		show_trace_task(p);
	}
}

char * render_sigset_t(sigset_t *set, char *buffer)
{
	int i = _NSIG, x;
	do {
		i -= 4, x = 0;
		if (sigismember(set, i+1)) x |= 1;
		if (sigismember(set, i+2)) x |= 2;
		if (sigismember(set, i+3)) x |= 4;
		if (sigismember(set, i+4)) x |= 8;
		*buffer++ = (x < 10 ? '0' : 'a' - 10) + x;
	} while (i >= 4);
	*buffer = 0;
	return buffer;
}

void show_state(void)
{
	struct task_struct *p;

#if (BITS_PER_LONG == 32)
	printk("\n"
	       "                         free                        sibling\n");
	printk("  task             PC    stack   pid father child younger older\n");
#else
	printk("\n"
	       "                                 free                        sibling\n");
	printk("  task                 PC        stack   pid father child younger older\n");
#endif
	read_lock(&tasklist_lock);
	for_each_task(p) {
		/*
		 * reset the NMI-timeout, listing all files on a slow
		 * console might take alot of time:
		 */
		touch_nmi_watchdog();
		show_task(p);
	}
	read_unlock(&tasklist_lock);
}

/**
 * reparent_to_init() - Reparent the calling kernel thread to the init task.
 *
 * If a kernel thread is launched as a result of a system call, or if
 * it ever exits, it should generally reparent itself to init so that
 * it is correctly cleaned up on exit.
 *
 * The various task state such as scheduling policy and priority may have
 * been inherited fro a user process, so we reset them to sane values here.
 *
 * NOTE that reparent_to_init() gives the caller full capabilities.
 */
void reparent_to_init(void)
{
	//printk(KERN_INFO "reparent_to_init (%d)\n", current->pid);
	write_lock_irq(&tasklist_lock);

	/* Reparent to init */
	REMOVE_LINKS(current);
	current->p_pptr = child_reaper;
	current->p_opptr = child_reaper;
	SET_LINKS(current);

	/* Set the exit signal to SIGCHLD so we signal init on exit */
	current->exit_signal = SIGCHLD;

	current->ptrace = 0;
	if ((current->policy == SCHED_OTHER) && (current->nice < DEF_NICE))
		set_user_nice(current, DEF_NICE);
	/* cpus_allowed? */
	/* rt_priority? */
	/* signals? */
	current->cap_effective = CAP_INIT_EFF_SET;
	current->cap_inheritable = CAP_INIT_INH_SET;
	current->cap_permitted = CAP_FULL_SET;
	current->keep_capabilities = 0;
	memcpy(current->rlim, init_task.rlim, sizeof(*(current->rlim)));
	current->user = INIT_USER;

	write_unlock_irq(&tasklist_lock);
	//printk(KERN_INFO "~reparent_to_init (%d)\n", current->pid);
}

/*
 *	Put all the gunge required to become a kernel thread without
 *	attached user resources in one place where it belongs.
 */

void daemonize(void)
{
	struct fs_struct *fs;


	/*
	 * If we were started as result of loading a module, close all of the
	 * user space pages.  We don't need them, and if we didn't close them
	 * they would be locked into memory.
	 */
	exit_mm(current);

	current->session = 1;
	current->pgrp = 1;
	current->tty = NULL;

	/* Become as one with the init task */

	exit_fs(current);	/* current->fs->count--; */
	fs = init_task.fs;
	current->fs = fs;
	atomic_inc(&fs->count);
 	exit_files(current);
	current->files = init_task.files;
	atomic_inc(&current->files->count);
}

extern unsigned long wait_init_idle;

static inline void double_rq_lock(runqueue_t *rq1, runqueue_t *rq2)
{
	if (rq1 == rq2)
		spin_lock(&rq1->lock);
	else {
		if (rq1->cpu < rq2->cpu) {
			spin_lock(&rq1->lock);
			spin_lock(&rq2->lock);
		} else {
			spin_lock(&rq2->lock);
			spin_lock(&rq1->lock);
		}
	}
}

static inline void double_rq_unlock(runqueue_t *rq1, runqueue_t *rq2)
{
	spin_unlock(&rq1->lock);
	if (rq1 != rq2)
		spin_unlock(&rq2->lock);
}

void __init init_idle(void)
{
	//printk(KERN_INFO "init_idle (%d)\n", current->pid);
	runqueue_t *this_rq = this_rq(), *rq = current->p_mlfq->rq;
	unsigned long flags;

	__save_flags(flags);
	__cli();
	double_rq_lock(this_rq, rq);

	this_rq->curr = this_rq->idle = current;
	deactivate_task(current, rq);
	current->p_mlfq = NULL;
	current->priority = MAX_PRIO;
	current->state = TASK_RUNNING;
	clear_bit(smp_processor_id(), &wait_init_idle);
	double_rq_unlock(this_rq, rq);
	while (wait_init_idle) {
		cpu_relax();
		barrier();
	}
	current->need_resched = 1;
	__sti();
	//printk(KERN_INFO "~init_idle\n");
}

extern void init_timervecs (void);

void __init sched_init(void)
{
	//printk(KERN_INFO "sched_init\n");
	int k, nr, cpu=0;
	runqueue_t *rq = cpu_rq(0);
	mlfq_t *p_mlfq;

	rq->p_mlfq = rq->mlfq_instant + 0;
	spin_lock_init(&rq->lock);

	rq->cpu = 0;
	rq->nr_running = 0;
	rq->nr_switches = 0;

	p_mlfq = rq->p_mlfq;
	p_mlfq->rq = rq;
	p_mlfq->lock = &rq->lock;
	for (k = 0; k < MAX_PRIO; k++) {
		INIT_LIST_HEAD(p_mlfq->queue + k);
		p_mlfq->bitmap[k] = 1;
	}

	init_task.processor = cpu;

	rq = this_rq();
	rq->curr = current;
	rq->idle = NULL;
	wake_up_process(current);

	for(nr = 0; nr < PIDHASH_SZ; nr++)
		pidhash[nr] = NULL;

	init_timervecs();

	init_bh(TIMER_BH, timer_bh);
	init_bh(TQUEUE_BH, tqueue_bh);
	init_bh(IMMEDIATE_BH, immediate_bh);

	/*
	 * The boot idle thread does lazy MMU switching as well:
	 */
	atomic_inc(&init_mm.mm_count);
	enter_lazy_tlb(&init_mm, current, cpu);
	//printk(KERN_INFO "~sched_init\n");
}
