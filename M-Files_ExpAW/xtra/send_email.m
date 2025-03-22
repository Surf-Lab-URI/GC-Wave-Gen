function send_email(mail,password,server, obj,content)

% Apply prefs and props
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port', '587');
props.setProperty('mail.smtp.starttls.enable','true');
props.setProperty('mail.smtp.auth','true');
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server',server);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% Send email
sendmail('fabio.addona1@gmail.com',obj,content)