import logging
import re
from nomadcore.match_highlighter import MatchHighlighter, ANSI


LOGGER = logging.getLogger(__name__)


class Annotator(object):
    def __init__(self, annotateFilename=None,
                 formatNameWidth=15, formatSourceWidth=20):
        self.matchHighlighter = MatchHighlighter()
        if annotateFilename is None:
            self.annotateFile = None
        else:
            LOGGER.info("writing annotated input to " + annotateFilename)
            self.annotateFile = open(annotateFilename, 'w')
            self.annotateFile.write(
                "# this file contains ANSI colors; use 'less -R' to view it\n")
        self._formatNameWidth = formatNameWidth
        self._formatSourceWidth = formatSourceWidth
        self._update_annotate_format()

    def _update_annotate_format(self):
        self._annotate_format = '%%%ds:%%04d %%%ds %%9s|%%s' % (
            self._formatSourceWidth, self._formatNameWidth)

    def set_formatSourceWidth(self, formatSourceWidth):
        self._formatSourceWidth = formatSourceWidth
        self._update_annotate_format()

    def set_formatNameWidth(self, formatNameWidth):
        self._formatNameWidth = formatNameWidth
        self._update_annotate_format()

    def annotate(self, minfo):
        # don't write to file unless it is there
        if not self.annotateFile:
            return 1

        # setup match label
        matchlabel = ''
        if minfo['coverageIgnore'] == 1:
            matchlabel = 'l_ign'
        elif minfo['coverageIgnore']:
            matchlabel = 'g_ign'
        elif minfo['match']:
            if minfo['matcher_does_nothing']:
                matchlabel='n_'
            if minfo['match']<2:
                matchlabel += 'p_'
            matchlabel += minfo['which_re']
        else:
            matchlabel = 'no'

        # check if there is pre-highlighted data in minfo
        highlighted = minfo.get('highlighted', None)
        if highlighted is None:
            # highlight line
            if minfo['coverageIgnore']:
                highlighted = self.matchHighlighter.highlight_minfo(
                    minfo, linecolor=ANSI.FG_BLUE)
            elif minfo['match']:
                if minfo['matcher_does_nothing']:
                    highlighted = self.matchHighlighter.highlight_minfo(
                        minfo, linecolor=ANSI.FG_MAGENTA)
                else:
                    highlighted = self.matchHighlighter.highlight_minfo(minfo)
            else:
                highlighted = minfo['fInLine']
        # shorted matcher- and source file names
        name = minfo['matcherName'][-self._formatNameWidth:] if minfo['matcherName'] else 'UNNAMED'
        defFile = minfo['defFile'][-self._formatSourceWidth:]
        self.annotateFile.write(self._annotate_format % (
            defFile, minfo['defLine'], name,
            matchlabel, highlighted,
        ))
