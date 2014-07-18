#!/bin/python

#ui mail

email_addresses = {'crc_admins':['bposwald@gmail.com','bcheldelin@uidaho.edu'],'all':["msettles@uidaho.edu,dnew@uidaho.edu,bcheldelin@uidaho.edu,boswald@uidaho.edu"],'admins':["boswald@uidaho.edu,shunter@gmail.com"],'test':["boswald@uidaho.edu"]}
smtp_server = 'm.outlook.com'
smtp_port = 587


# Now construct the message
import smtplib, email
from email.mime.text import MIMEText
import netrc
import argparse

 
parser = argparse.ArgumentParser(description='Sends a simple email to predefined lists of email adresses')
parser.add_argument('-s','--subject', help='subject of the email',required=False, default="subject")
parser.add_argument('-m','--message',help='message of the email', required=False, default="an empty message")
parser.add_argument('-t','--to',help='one of: all, admins, crc_admins. Or a list of addresses separated by semicolons', required=False, default="test")
args = parser.parse_args()

#-----------------------------------------------------------------------------
# this function requires a .netrc file in the user's home directory of format:
# machine smtpserver.example.com
#      login  user.name@example.com
#      password somePassword

credentials = False

try:
	credentials = netrc.netrc()
	username, account, password = credentials.authenticators(smtp_server)
except IOError as error:
	print("Warning: no .netrc file found, emails aren't going to be sent - but will print to screen")

#-----------------------------------------------------------------------------

def send(recipients,subject,message): 
	if credentials:
		msg = MIMEText(message)

		msg['Subject'] = subject
		msg['From'] = username
		msg['To'] = ";".join(recipients)

		# Now send the message
		mailer = smtplib.SMTP(smtp_server, smtp_port)
		mailer.ehlo()
		mailer.starttls()
		mailer.ehlo()
		mailer.login(username,password)
		mailer.sendmail(username, recipients, msg.as_string())
		mailer.quit()
	else:
		print("|------------------email---------------------\n| "+subject+"\n| ----\n| "+message+"\n|--------------------------------------------\n")

if args.to in email_addresses:
	send(email_addresses[args.to],args.subject,args.message)
else:
	send(args.to.split(';'),args.subject,args.message)
