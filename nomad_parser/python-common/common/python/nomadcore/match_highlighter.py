from builtins import range
from builtins import object
import logging
from   operator import itemgetter
import re


LOGGER = logging.getLogger("nomadcore.match_highlighter")


class ANSI(object):
    RESET='\033[0m'
    RESET_COLOR='\033[39;49m'
    BEGIN_INVERT='\033[7m'
    END_INVERT='\033[27m'
    FG_BLACK='\033[30m'
    FG_RED='\033[31m'
    FG_GREEN='\033[32m'
    FG_YELLOW='\033[33m'
    FG_BLUE='\033[34m'
    FG_MAGENTA='\033[35m'
    FG_CYAN='\033[36m'
    FG_WHITE='\033[37m'
    BG_BLACK='\033[40m'
    BG_RED='\033[41m'
    BG_GREEN='\033[42m'
    BG_YELLOW='\033[43m'
    BG_BLUE='\033[44m'
    BG_MAGENTA='\033[45m'
    BG_CYAN='\033[46m'
    BG_WHITE='\033[47m'
    FG_BRIGHT_BLACK='\033[30;1m'
    FG_BRIGHT_RED='\033[31;1m'
    FG_BRIGHT_GREEN='\033[32;1m'
    FG_BRIGHT_YELLOW='\033[33;1m'
    FG_BRIGHT_BLUE='\033[34;1m'
    FG_BRIGHT_MAGENTA='\033[35;1m'
    FG_BRIGHT_CYAN='\033[36;1m'
    FG_BRIGHT_WHITE='\033[37;1m'
    BG_BRIGHT_BLACK='\033[40;1m'
    BG_BRIGHT_RED='\033[41;1m'
    BG_BRIGHT_GREEN='\033[42;1m'
    BG_BRIGHT_YELLOW='\033[43;1m'
    BG_BRIGHT_BLUE='\033[44;1m'
    BG_BRIGHT_MAGENTA='\033[45;1m'
    BG_BRIGHT_CYAN='\033[46;1m'
    BG_BRIGHT_WHITE='\033[47;1m'


RE_finalnewlines=re.compile(r"^(.*?)(\n*)$")

class MatchHighlighter(object):
    """ANSI-Color highlighting for capturing groups in regex matches"""
    def __init__(self, multiline=False):
        self.multiline = multiline

    def highlight_minfo(self, minfo, multiline=None, linecolor=''):
        """
        :rtype: str
        :returns:  highlighed string containing ANSI escapes
        """
        if multiline is None:
           multiline=self.multiline

        line = minfo['fInLine']
        trailing = ''
        if not multiline:
            m = RE_finalnewlines.match(line)
            line = m.group(1)
            trailing = m.group(2)

        # list with the highlight switching events
        # tuples: (pos, event, ansi)
        #   event:
        #     matchstart=-1
        #          close= 0
        #           open= 1
        #       matchend= 2
        #    ansi: escape code
        ansiSwitch = list()
        # mark whole area of match by fg/bg inversion
        ansiSwitch.append([minfo['span'][0][0][0],-1,ANSI.BEGIN_INVERT])
        ansiSwitch.append([minfo['span'][0][0][1], 2,ANSI.END_INVERT])
        # capture groups
        nspan=len(minfo['span']) # we covered 0 (overall capture) before...
        for groupi in range(1,nspan):
            for s,e in minfo['span'][groupi]:
                ansiSwitch.append([s,1,ANSI.FG_RED])
                ansiSwitch.append([e,0,ANSI.RESET_COLOR])
        # sort by position, then by event
        ansiSwitch=sorted(ansiSwitch, key=itemgetter(0,1))
        if ansiSwitch[-1][0]>len(minfo['fInLine']):
            raise Exception('line_exceed',"match exceeds fInLine")
        # line may not contain trailing newline characters
        for a in ansiSwitch:
            if a[0] > len(line):
                a[0] = len(line)
        # append line remainder and always terminate by ANSI reset
        ansiSwitch.append([len(line),3,ANSI.RESET])
        lastPos=0
        highlighted=linecolor
        for i in range(len(ansiSwitch)):
            if ansiSwitch[i][0]>lastPos:
                highlighted += line[lastPos:ansiSwitch[i][0]]
                lastPos = ansiSwitch[i][0]
            highlighted += ansiSwitch[i][2]
        # re-append trailing newlines if we stripped them
        highlighted += trailing
        return highlighted

    def highlight(self, match, fInLine, multiline=None, linecolor=''):
        """
        :rtype: str
        :returns:  highlighed string containing ANSI escapes
        """
        # subset of annotator.match_info for backwards compatibility
        minfo = {
            'fInLine': fInLine,
            'span': [[match.span()]],
        }
        ngroups=len(match.groups())
        for groupi in range(1,ngroups+1):
            if match.group(groupi) is None:
                minfo['span'].append([])
            else:
                minfo['span'].append([match.span(groupi)])
        return self.highlight_minfo(minfo, multiline, linecolor)
