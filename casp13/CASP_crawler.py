#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Jun-3, latest modified on 2018-Jun-28


import urllib.request
from bs4 import BeautifulSoup
import datetime
import requests
import argparse
import time
import sys


# Download HTML
def get_HTML_content(url):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/55.0.2883.87 Safari/537.36'}  # Add header
    req = urllib.request.Request(url, headers=headers)
    html = urllib.request.urlopen(req)
    content = html.read().decode('utf-8')  # Transcode to utf-8
    html.close()  # Remember to close the website, otherwise you will get error
    return content


def get_namelist(content):
    namelist = []
    print ("The list of regular protein in CASP13 at current time is ")
    soup = BeautifulSoup(content, 'html.parser')
    tr = soup.find_all('tr', {'class': 'datarow'})
    for elements in tr:
        print(elements.find_all('td')[1].contents[0].get_text())
        # print(elements.find('a').get_text()) alternative way to get target name list
        namelist += elements.find_all('td')[1].contents[0].get_text()
    return tr


def check_today(tr, today):
    linklist = []
    for elements in tr:
        if elements.find_all('td')[5].contents[0] == today:
            name = elements.find_all('td')[1].contents[0].get_text()
            pesudolink = elements.find_all('td')[1].contents[0].get('href')
            link = "http://predictioncenter.org/casp13/" + pesudolink
            print ("Today's protein is %s, its link is %s" % (name, link))
            linklist.append(link)
    return linklist


def get_fasta(linklist, email):
    for link in linklist:
        soup = BeautifulSoup(get_HTML_content(link), 'html.parser')
        td = soup.find_all('td', {'colspan': '2'})
        name = td[0].contents[0].get_text().split(': ')[1]
        fasta = td[1].find('textarea', {'class': 'input_pdb'}).contents[0]
        with open('%s.fasta' % name, 'w') as fw:
            fw.write(fasta)
            fw.close()
        print(fasta)
        oddeven = name[1:].split('s')
        print(oddeven)
        if len(oddeven) == 2:
            pound = oddeven[1]
        else:
            pound = oddeven[0]
        if len(email) > 1:
            submit(name, email[int(pound) % len(email)], fasta)
        else:
            submit(name, email, fasta)


def submit(name, email, fasta):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/45.0.2454.85 Safari/537.36',
    }

    r = requests.post('http://raptorx.uchicago.edu/ContactMap/',
                      files={'jobname': (None, name), 'email': (None, email),
                             'sequences': (None, fasta),
                             'fileseq': (None, ''), 'predict_sub': (None, 'Submit')}, headers=headers)

    soup = BeautifulSoup(r.text, 'html.parser')
    joblink = soup.find('h5').contents[1].get('href')
    print("The submitted job link is %s" % joblink)
    sys.stdout.flush()


def main():
    parser = argparse.ArgumentParser(
        description="This script submits today's CASP regular targets to RaptorX automatically")
    parser.add_argument("-d", "--date", help="date in yyyy-mm-dd format, default is today", type=str,
                        default=datetime.date.today())
    parser.add_argument("-e", "--email", nargs='+', help="1 or more email addresses to receive RaptorX", type=str,
                        default='1964569056@qq.com')
    args = parser.parse_args()
    today = args.date
    email = args.email

    today = str(today)[0:10]
    print("The day you choose is %s" % today)
    url = 'http://predictioncenter.org/casp13/targetlist.cgi?view=regular'
    html = get_HTML_content(url)
    protein = get_namelist(html)
    print("\n")
    linklist = []
    print("Time is %s" % datetime.datetime.now().strftime('%H:%M:%S'))
    linklist = check_today(protein, today)
    sys.stdout.flush()  # this command is used to print text immediately, not after program finished
    while linklist == []:
        time.sleep(15)
        html = get_HTML_content(url)
        protein = get_namelist(html)
        print("\n")
        linklist = check_today(protein, today)
        print("Time is %s" % datetime.datetime.now().strftime('%H:%M:%S'))
        sys.stdout.flush()

    get_fasta(linklist, email)


if __name__ == '__main__':
    main()
