from fireworks import LaunchPad

def create_launchpad(username, password, server = "serenity", lpadname = None):
    """
    Creates the fireworks launchpad on specific preset servers.

    Args:
        username (str) : username for the mongodb database
        password (str) : password for the mongodb database
        server   (str) : server name: "serinity" (default) or "atlas"
        lpadname (str) : name of the fireworks internal database. If not given,
                         the name is inferred.

    Returns:
        fireworks object : Launchpad for internal fireworks use.
    """
    if server == "atlas":
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
