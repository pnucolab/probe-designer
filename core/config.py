import os
import uuid
task_id = str(uuid.uuid4())
os.makedirs("./files", exist_ok=True)
target_tarnscript =  f"./files/{task_id}_output.fa"
header_file = f"./files/{task_id}_target_transcript.txt"