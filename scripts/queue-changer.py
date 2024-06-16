import subprocess
import time
import sys
import datetime

def get_job_status(job_id: int):
    try:
        result = subprocess.run(["squeue", "--job", str(job_id)], capture_output=True, text=True, check=True)
        
        lines = result.stdout.strip().split('\n')
        if len(lines) > 1:
            # Skip the header line
            status_line = lines[1]
            job_status = status_line.split()[4]
            return job_status
        else:
            return "NOT_FOUND"
    except subprocess.CalledProcessError:
        return "CalledProcessError"
      
def change_queue_submission(job_id: int, queue_name: str):
    try:
        subprocess.run(["scontrol", "update", f"jobid={job_id}", f"partition={queue_name}"])
    except subprocess.CalledProcessError:
        return "CalledProcessError"


def monitor_job(job_id: int, queue_name: str):
    while True:
        job_status = get_job_status(job_id)
        print(f"Job {job_id} status: {job_status}")
        
        if job_status == "CG":
            print(f'Revocation: {datetime.datetime.now()}')
            change_queue_submission(job_id, queue_name)
        elif job_status in ["NOT_FOUND", "CalledProcessError"]:
            break

        time.sleep(30)

if __name__ == "__main__":
    job_id = int(sys.argv[1])
    queue_name = str(sys.argv[2])
    monitor_job(job_id, queue_name)