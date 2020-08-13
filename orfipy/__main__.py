#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:57:38 2020

@author: usingh
"""
import argparse

def main():
    parser = argparse.ArgumentParser(description="""
     
    description='orfipy: extract Open Reading Frames',
    usage='''orfipy [<options>] <infasta>
    """)
    parser.add_argument("name", help="Name of Employee")
    parser.add_argument("title", help="Job Title of Employee")
    parser.add_argument("--address", help="Address of Employee")
    parser.add_argument('infile', help='Fasta containing sequences',action="store")
    args = parser.parse_args()
    

      
    
    
    
if __name__ == '__main__':
    main()