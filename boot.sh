export CONFIG_OBJ='config.ProductionConfig'
adduser --disabled-password --gecos '' r
cd /creeds/
mod_wsgi-express setup-server wsgi.py --port=80 --user r --group r --server-root=/etc/creeds --socket-timeout=600
chmod a+x /etc/creeds/handler.wsgi
/etc/creeds/apachectl start
tail -f /etc/creeds/error_log
