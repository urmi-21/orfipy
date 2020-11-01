#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 22:11:51 2020

@author: usingh
"""

import sys
import datetime as dt
from colorama import Fore,Style


#color print
def orfipy_print(color,*args,stderr=False,**kwargs):
    if stderr:
        print(color,file=sys.stderr,end="",**kwargs)
        print(*args,file=sys.stderr,end="",**kwargs)
        print(Style.RESET_ALL,file=sys.stderr,**kwargs)
    else:
        print(color,file=sys.stdout,end="",**kwargs)
        print(*args,file=sys.stdout,end="",**kwargs)
        print(Style.RESET_ALL,file=sys.stdout,**kwargs)
def print_notification(*args):
    """Print message to stderr
    """
    orfipy_print(Fore.LIGHTYELLOW_EX,*args,stderr=True)
def print_message(*args):
    """Print message to stderr
    """
    orfipy_print(Fore.LIGHTCYAN_EX,*args,stderr=True)
def print_error(*args):
    """Print message to stderr
    """
    orfipy_print(Fore.LIGHTRED_EX,*args,stderr=True)
    
def print_success(*args):
    """Print message to stderr
    """
    orfipy_print(Fore.LIGHTGREEN_EX,*args,stderr=True)
    
def get_time_stamp():
    """Function to return current timestamp.
    Parameters
    ----------
    shorten: bool
        return short version without space, dash and colons
    :return: timestamp as string
    :rtype: string
    """
    timestamp=str(dt.datetime.now()).replace(" ","-").replace("-","_").replace(" ","").replace(":","_")
    return timestamp