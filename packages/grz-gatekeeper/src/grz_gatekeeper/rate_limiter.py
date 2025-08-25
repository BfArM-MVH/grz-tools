import threading
import time

from fastapi import HTTPException, Request, status

# FIXME: ideally we'd use some existing solution such as slowapi instead
RATE_LIMIT_CACHE = {}
RATE_LIMIT_LOCK = threading.Lock()


class RateLimiter:
    def __init__(self, times: int, seconds: int):
        self.times = times
        self.seconds = seconds

    def __call__(self, request: Request):
        ip = request.client.host if request.client else "unknown"
        with RATE_LIMIT_LOCK:
            current_time = time.time()

            if ip in RATE_LIMIT_CACHE:
                last_time, count = RATE_LIMIT_CACHE[ip]
                if current_time - last_time > self.seconds:
                    RATE_LIMIT_CACHE[ip] = (current_time, 1)
                else:
                    if count >= self.times:
                        raise HTTPException(
                            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                            detail="Too Many Requests",
                        )
                    RATE_LIMIT_CACHE[ip] = (last_time, count + 1)
            else:
                RATE_LIMIT_CACHE[ip] = (current_time, 1)
