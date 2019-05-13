#!/usr/bin/env python3
#
# Interactive display of ancestry using xdot from Jose Fonseca.
# https://github.com/jrfonseca/xdot.py
#

import sys
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
import xdot

input_file = str(sys.argv[1])
dotgraph = open(input_file).read()


dotcode = bytes(dotgraph, "utf-8")

class MyDotWindow(xdot.DotWindow):

    def __init__(self):
        xdot.DotWindow.__init__(self)
        self.dotwidget.connect('clicked', self.on_url_clicked)

    def on_url_clicked(self, widget, url, event):
        dialog = Gtk.MessageDialog(
            parent=self,
            buttons=Gtk.ButtonsType.OK,
            message_format="%s" % url)
        dialog.connect('response', lambda dialog, response: dialog.destroy())
        dialog.run()
        return True

def main():
    window = MyDotWindow()
    window.set_dotcode(dotcode)
    window.connect('delete-event', Gtk.main_quit)
    Gtk.main()


if __name__ == '__main__':
    main()

