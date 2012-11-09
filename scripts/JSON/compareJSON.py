#!/usr/bin/env python

import sys
import optparse
from FWCore.PythonUtilities.LumiList import LumiList


if __name__ == '__main__':
    
    parser = optparse.OptionParser ("Usage: %prog --command [--options] alpha.json beta.json [output.json]")
    # required parameters
    cmdGroup = optparse.OptionGroup (parser, "Command Options ")
    cmdGroup.add_option ('--and', dest='command', action='store_const',
                         const='and', 
                         help = '"and" (i.e., take intersection) of two files')
    cmdGroup.add_option ('--or', dest='command', action='store_const',
                         const='or',
                         help = '"or" (i.e., take union) of two files')
    cmdGroup.add_option ('--sub', dest='command', action='store_const',
                         const='sub',
                         help = '"subtraction" (i.e., lumi sections in alpha not in beta) of two files')
    cmdGroup.add_option ('--diff', dest='command', action='store_const',
                         const='diff',
                         help = 'All differences (i.e., alpha - beta AND beta - alpha) of two files. Output will only be to screen (not proper JSON format).')
    parser.add_option_group (cmdGroup)
    (options, args) = parser.parse_args()
    if len (args) < 2:
        raise RuntimeError, "At least two input filenames must be provided."
    if not options.command:
        raise RuntimeError, "Exactly one command option must be specified"

    LumiListList=[]
    for arg in args: LumiListList.append(LumiList(arg))

    ##################
    ## Diff Command ##
    ##################
    if options.command == 'diff':
        if len (args) >= 3:
            raise RuntimeError, "Can not output to file with '--diff' option.  The output is not standard JSON."
        firstOnly  = LumiListList[0] - LumiListList[1]
        secondOnly = LumiListList[1] - LumiListList[0]
        if not firstOnly and not secondOnly:
            print "Files '%s' and '%s' are the same." % (args[0], args[1])
            sys.exit()
        print "'%s'-only lumis:" % args[0]
        if firstOnly:
            print firstOnly
        else:
            print "None"
        print "\n'%s'-only lumis:" % args[1]
        if secondOnly:
            print secondOnly
        else:
            print "None"
        sys.exit()
    
    ########################
    ## All other commands ##
    ########################
    if options.command == 'and':
        outputList = LumiListList[0]
        for lumi in LumiListList[1:]:
            outputList = outputList & lumi

    if options.command == 'or':
        outputList = LumiListList[0]
        for lumi in LumiListList[1:]:
            outputList = outputList | lumi

    if options.command == 'sub':
        outputList = LumiListList[0]
        for lumi in LumiListList[1:]:
            outputList = outputList - lumi

    print outputList
