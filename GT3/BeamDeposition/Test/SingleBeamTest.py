#!/usr/bin/python

from TestBase.GT3TestBase import GT3TestBase
from GT3.gt3 import gt3

if __name__=="__main__":

    try:
        test = GT3TestBase()
        shot = gt3(preparedInput=test)
        shot.run_NBI()
        print("Test successful")
    except Exception as e:
        print("Test failed")
        raise e



