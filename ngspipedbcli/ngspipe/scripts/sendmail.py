import os
import argparse
import sys
import time
import datetime
from smtplib import SMTP
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.header import Header
 
parser = argparse.ArgumentParser(description="Sending emails with attachments to users.") 
parser.add_argument('-r', '--receiver', type=str, required=True, help='multiple email address with dogma inside', default = '')
parser.add_argument('-s', '--sender', type=str, default='ngspipedb@sina.com', required=False, help='not nessary')
parser.add_argument('-p', '--password', type=str, default='8554ee2faf7409fe', required=False, help='not nessary')
parser.add_argument('-smtp', '--smtpserver', type=str, default='smtp.sina.com.cn', required=False, help='not nessary')
parser.add_argument('-t', '--status', type=str, default='success', required=True, help='success or error')
parser.add_argument('-d', '--logdir', type=str, required=False, help='not nessary')
parser.add_argument('-l', '--logfile', type=str, required=False, help='not nessary')


def new_report(logdir):
    logfiles = os.listdir(logdir)                                    #列出目录的下所有文件和文件夹保存到lists
    logfiles.sort(key=lambda fn:os.path.getmtime(os.path.join(logdir, fn)))#按时间排序
    file_new = os.path.join(logdir,logfiles[-1])                     #获取最新的文件保存到file_new
    return file_new

if __name__ == "__main__":

    #Email content template
    SendHtml = """
    <html>
        <head>
        </head>
        <body background-color='red'>
        <p>{}</p>
        <p>Dear user, </p>
        <p>Please find attached log file.</p>
        </p>
            if you like this project?<br>
            <a href="https://github.com/xuanblo/NGSPipeDb">NGSPipeDb</a> 
            Please star it on Github.
        </p>
        </body>
    </html>
    """.format(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))

    args = parser.parse_args()

    # attach file
    if args.logfile:
        attachmentFilePath = args.logfile
    if args.logdir:
        attachmentFilePath = new_report(args.logdir)

    sendfile = open(attachmentFilePath,"rb").read()
    addfile = MIMEText(sendfile,"base64","utf-8")
    addfile["Content-Type"] = "application/octet-stream"
    addfile["Content-Disposition"] = "attachment; filename={}".format(os.path.basename(attachmentFilePath))
    
    # main information
    sender = args.sender
    password = args.password
    receiver = args.receiver

    if '@' not in receiver or receiver == "":
        sys.stderr.write('No email address is provided.\n')
        sys.exit(0)

    subject = "[NGSPipeDb] - Your analysis: {}".format(args.status.upper())
    msg = MIMEMultipart("related")
    msg["Subject"]=Header(subject,"utf-8")
    msg['From'] = Header(sender) #发件人
    msg['to'] = Header(receiver)  #收件人（主送）
    msgText = MIMEText(SendHtml,"html","utf-8")

    msg.attach(msgText) #添加邮件正文到邮件中
    msg.attach(addfile) #添加附件到邮件中

    ret = True
    try:
        # Create secure connection with server and send email
        smtp = SMTP()
        smtp.connect(args.smtpserver)
        smtp.login(sender, password)
        smtp.sendmail(from_addr=sender, 
                        to_addrs=receiver.split(','), 
                        msg=msg.as_string())
        smtp.quit() 
    except Exception:
        ret = False

    if ret:
        sys.stderr.write("Mail was successfully sended.\n")
    else:
        sys.stderr.write("Mail was failed sended.\n")


