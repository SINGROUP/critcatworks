from fabric import Connection

FW_DIRECTORY = ".fireworks"
REMOTE_HOST = 'jagermar@taito.csc.fi'
PYTHON_INSTRUCTIONS = "launch_instructions.py"

def remote_launch():
    c = Connection(REMOTE_HOST)
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