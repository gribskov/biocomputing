#!/usr/bin/env python3
# This script performs a GET and a POST request through HTTP/HTTPS using
# builtin python3 moudules. There is also a class to encode files for upload!

import urllib.request
import http.client
import mimetypes
import codecs
import uuid
import binascii
import io
import os
import sys

REMOTE_PROTOCOL = 'http' # or https
REMOTE_HOST = 'httpbin.org'
REMOTE_PORT = '80' # 443 for https
USERAGENT = 'MySuperUserAgent'

class MultipartFormdataEncoder(object):
    def __init__(self):
        self.boundary = uuid.uuid4().hex
        self.content_type = 'multipart/form-data; boundary={}'.format(self.boundary)

    @classmethod
    def u(cls, s):
        if sys.hexversion < 0x03000000 and isinstance(s, str):
            s = s.decode('utf-8')
        if sys.hexversion >= 0x03000000 and isinstance(s, bytes):
            s = s.decode('utf-8')
        return s

    def iter(self, fields, files):
        """
        fields is a sequence of (name, value) elements for regular form fields.
        files is a sequence of (name, filename, file-type) elements for data to be uploaded as files
        Yield body's chunk as bytes
        """
        encoder = codecs.getencoder('utf-8')
        for (key, value) in fields:
            key = self.u(key)
            yield encoder('--{}\r\n'.format(self.boundary))
            yield encoder(self.u('Content-Disposition: form-data; name="{}"\r\n').format(key))
            yield encoder('\r\n')
            if isinstance(value, int) or isinstance(value, float):
                value = str(value)
            yield encoder(self.u(value))
            yield encoder('\r\n')
        for (key, filename, fpath) in files:
            key = self.u(key)
            filename = self.u(filename)
            yield encoder('--{}\r\n'.format(self.boundary))
            yield encoder(self.u('Content-Disposition: form-data; name="{}"; filename="{}"\r\n').format(key, filename))
            yield encoder('Content-Type: {}\r\n'.format(mimetypes.guess_type(filename)[0] or 'application/octet-stream'))
            yield encoder('\r\n')
            with open(fpath,'rb') as fd:
                buff = fd.read()
                yield (buff, len(buff))
            yield encoder('\r\n')
        yield encoder('--{}--\r\n'.format(self.boundary))

    def encode(self, fields, files):
        body = io.BytesIO()
        for chunk, chunk_len in self.iter(fields, files):
            body.write(chunk)
        return self.content_type, body.getvalue()


if len(sys.argv) < 2:
    sys.stderr.write('Please pass a file as argument\n')
    sys.exit(1)

if not os.path.isfile(sys.argv[1]):
    sys.stderr.write('%s is not a valid file\n' % sys.argv[1])
    sys.exit(1)


# SIMPLE HTTP GET
req = urllib.request.Request(REMOTE_PROTOCOL+'://'+REMOTE_HOST+':'+REMOTE_PORT+'/get',headers={'User-Agent':USERAGENT})
with urllib.request.urlopen(req,timeout=10) as f:
    buf = f.read(300).decode('utf-8')
    print("GET RESPONSE")
    print(buf)

# POST FIELDS AND FILES
fields = [('param1','paramvalue1'),('param2','paramversion2'),('param3','paramvalue3')]
files = [('fieldname',os.path.basename(sys.argv[1]),sys.argv[1])]
content_type, body = MultipartFormdataEncoder().encode(fields, files)
h = None
if REMOTE_PROTOCOL == 'http':
    h = http.client.HTTPConnection(REMOTE_HOST,REMOTE_PORT)
else:
    h = http.client.HTTPSConnection(REMOTE_URL,REMOTE_PORT)
headers = {
  'User-Agent': USERAGENT,
  'Content-Type': content_type
}
h.request('POST', '/post', body, headers)
res = h.getresponse()
if not res.status in (200,201):
    sys.stderr.write(res.read().decode('utf-8')+'\n')
    sys.exit(1)
print("POST RESPONSE")
print(res.read().decode('utf-8'))