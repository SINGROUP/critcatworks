from fireworks import LaunchPad

def create_launchpad(username, password, server = "serenity"):
    if server == "atlas":
        if username == "mjcritcat":
            name = "fireworks"
        else:
            name = username[:2] + "fireworks"

        lp = LaunchPad(host = "austerity-shard-00-01-hgeov.mongodb.net:27017",
            port = 27017,
            name = name,
            username = username,
            password = password,
            logdir =  ".",
            strm_lvl = "INFO",
            ssl =  True,
            authsource = "admin")
    elif server == "serenity":
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
