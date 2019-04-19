from fireworks import LaunchPad

def create_launchpad(username, password, server = "serenity", lpadname = None):
    if server == "atlas":
        name = username[:2] + "fireworks"

        lp = LaunchPad(host = "austerity-shard-00-00-hgeov.mongodb.net:27017",
            port = 27017,
            name = name,
            username = username,
            password = password,
            logdir =  ".",
            strm_lvl = "INFO",
            ssl =  True,
            authsource = "admin")
    elif server == "serenity":
        if lpadname:
            name = lpadname
        else:
            name = username[:2] + "fireworks"
        lp = LaunchPad(host = "nanolayers.dyndns.org:27017",
        port = 27017,
        name = name,
        username = username,
        password = password,
        logdir =  ".",
        strm_lvl = "INFO",
        #ssl =  True,
        authsource = name)
    else:
        lp = LaunchPad()

    return lp
