from fireworks import LaunchPad

def create_launchpad():
    lp = LaunchPad(host = "austerity-shard-00-00-hgeov.mongodb.net:27017",
        port = 27017,
        name = "fireworks",
        username = "mjcritcat",
        password = "heterogeniuscatalysis",
        logdir =  ".",
        strm_lvl = "INFO",
        ssl =  True,
        authsource = "admin")

    return lp