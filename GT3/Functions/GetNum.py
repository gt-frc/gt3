#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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
                return float(input(prompt))
            except:
                print("Invalid input \n")
    elif valType=='i':
        while True:
            try:
                x = input(prompt)
                if isinstance(x, (int, long)):
                    return(int(x))
                else:
                    print "Invalid input \n"
                    continue
            except:
                print "Invalid input \n"