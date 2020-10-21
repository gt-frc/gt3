#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from deprecation import deprecated as Deprecated

@Deprecated(deprecated_in="0.0.3", details="Interactive mode is to be removed in 0.0.4", removed_in="0.0.4")
def getNum(prompt, valType):
    """

    Function to obtain valid inputs from keyboard

    :param prompt: Text prompt to user
    :param valType: Type of value expected from user
    :return: Your mom
    """

    if valType=='f':
        while True:
            try:
                return float(eval(input(prompt)))
            except:
                print("Invalid input \n")
    elif valType=='i':
        while True:
            try:
                x = eval(input(prompt))
                if isinstance(x, int):
                    return(int(x))
                else:
                    print("Invalid input \n")
                    continue
            except:
                print("Invalid input \n")