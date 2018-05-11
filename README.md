# linux-2.4.27-CS518
Implement multilevel feedback queues for scheduling in kernel 2.4.27

##Scheduler:
1. Currently there are 256 priority queues, with priority range from [0, 256). Time slices for the queues vary from 10 ms to 300 ms.
2. If a task yields before expiration of time-slice then it is simple placed at the tail of same queue. We check this in schedule(), by checking the value of remaining counter value, if it is not zero, then it means it yielded, before expiration of its time slice.
3. If the time-slice expires, then task is put at the tail of the next priority queue. We make this check in function "handle_tick_process()", which we call from the "do_timer()"->"update_process_times()", which is called every HZ ms. 
4. If a task is waiting for some I/O, then we shift it to a higher priority queue. We make this change when a task is put in wait queue and consequently removed from scheduler's queue in function "schedule()".

##Priority Inversion:
1. We implemented priority inversion using a technique called priority inheritance.
2. So if a high priority process needs to get a lock held by a low priority process, then we bump the priority of low-priority process to match that of high-priority process. This is done so that the low priority task can release the resource as soon as possible.
3. After the lock is released we check; if the priority was bumped, if so we change it back to original priority.
4. We implemented this only for the case of mutexes and not semaphores. 

