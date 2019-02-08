from fabric import Connection
from fabric import task

FW_DIRECTORY = ".fireworks"
REMOTE_HOST = 'jagermar@taito.csc.fi'
#REMOTE_HOST = 'jagerm1@triton.aalto.fi'
PYTHON_INSTRUCTIONS = "execute_firework_jobs.py"

def remote_launch():
    #connect_kwargs={
    #    "key_filename": "/home/myuser/.ssh/id_rsa_taito",
    #     }
    #c = Connection(REMOTE_HOST, connect_kwargs)
    c = Connection(REMOTE_HOST, )
    result = c.run('echo "Hello, world!"')
    result = c.run('mkdir -p .fireworks')
    #result = c.run('mkdir -p .fireworks/fw_logs')

    # cd currently works only with chaining
    #result = c.run('cd .fireworks')

    # yaml not needed for instructions
    #c.put('austerity_launchpad.yaml', FW_DIRECTORY + "/" )
    
    # instructions given in python
    c.put(PYTHON_INSTRUCTIONS, FW_DIRECTORY + "/")

    result = c.run('python3 ' + FW_DIRECTORY + "/" + PYTHON_INSTRUCTIONS)


if __name__ == '__main__':
    remote_launch()
