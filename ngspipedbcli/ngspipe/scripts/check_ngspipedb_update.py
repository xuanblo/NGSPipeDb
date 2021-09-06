import sys
import os
import urllib.request
import markdown
from bs4 import BeautifulSoup

version = sys.argv[1]

change_log_in_github = "https://raw.githubusercontent.com/xuanblo/NGSPipeDb/master/mkdocs/markdowns/changelog.md"
change_log_in_github = "http://www.liu-lab.com/ngspipedb/changelog.md"
#change_log_in_github = "/tmp/changelog.md"

def parse_with_wget(change_log_in_markdown, version):
    os.system(r"wget {} -O /tmp/changelog.md".format(change_log_in_markdown))

    with open("/tmp/changelog.md", 'r') as f:
        content = f.read()
        #print(content)
        
        html = markdown.markdown(content)
        #print(html)
        soup = BeautifulSoup(html, 'html.parser')
        #print(soup)
        targets = soup.find_all('h2')
        versions = [i.text[1:6] for i in targets]
        #print(versions)
        current_version_index = versions.index(version)
        
        if current_version_index:
            print("[Your pipeline need update] the lastest version is {}".format(versions[0]))
            print("[Your pipeline need update] your current version is {}".format(version))
            print("[Your pipeline need update] Please see https://xuanblo.github.io/NGSPipeDb/NGSPipe-RNA-seq/#NGSPipeDbSource")
            for target in targets[0:current_version_index]:
                print("[Your pipeline need update] ## {}".format(target.text))
                for sib in target.find_next_siblings():
                    if sib.name=="h2":
                        break
                    else:
                        update_text = sib.text.strip()
                        for line in update_text.split('\n'):
                            print("[Your pipeline need update] {}".format(line))
        else:
            print("[Your pipeline is the lastest version]")

def parse_with_urllib(change_log_in_markdown, version):
    with urllib.request.urlopen(change_log_in_markdown) as f:
        content = f.read().decode('utf-8')
        html = markdown.markdown(content)
        #print(html)
        soup = BeautifulSoup(html, 'html.parser')
        #print(soup)
        targets = soup.find_all('h2')
        versions = [i.text[1:6] for i in targets]
        #print(versions)
        current_version_index = versions.index(version)
        message = ""
        if current_version_index:
            message += "[Your pipeline need update] the lastest version is {}\n".format(versions[0])
            message += "[Your pipeline need update] your current version is {}\n".format(version)
            message += "[Your pipeline need update] Please see https://xuanblo.github.io/NGSPipeDb/NGSPipe-RNA-seq/#NGSPipeDbSource\n"
            for target in targets[0:current_version_index]:
                message += "[Your pipeline need update] ## {}\n".format(target.text)
                for sib in target.find_next_siblings():
                    if sib.name=="h2":
                        break
                    else:
                        update_text = sib.text.strip()
                        for line in update_text.split('\n'):
                            message += "[Your pipeline need update] {}\n".format(line)
        else:
            message += "[Your pipeline is the lastest version]\n"

    return message

#parse_with_wget(change_log_in_github, version)
#parse_with_urllib(change_log_in_github, version)

def ask_if_update_auto():
    return 0
    a = input('[oh-my-zsh] Would you like to update? [Y/n] y')
