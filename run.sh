#!/bin/bash
## For nginx-gunicorn-flask stack
# /etc/init.d/nginx start
# rm /etc/nginx/sites-enabled/default
# ln -s /my_application/creeds.conf /etc/nginx/sites-enabled/creeds

# nginx -t
# /etc/init.d/nginx restart

# cd /my_application && gunicorn wsgi:myapp \
# 	-b 127.0.0.1:8000 \
# 	-w 4 \
# 	--error-logfile=- \
# 	--access-logfile=- 

## For supervisor-nginx-gunicorn-flask stack
cp /my_application/supervisord.conf /etc/
cp /my_application/nginx.conf /etc/nginx/

service nginx stop

supervisord -c /etc/supervisord.conf -n
