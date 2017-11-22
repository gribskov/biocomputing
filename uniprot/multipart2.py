import http.client, mimetypes

def post_multipart(host, selector, fields, files):
    """
    Post fields and files to an http host as multipart/form-data.
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for data to be uploaded as files
    Return the server's response page.
    """
    content_type, body = encode_multipart_formdata(fields, files)

    # Choose between http and https connections
    if(selector.find('https') == 0):
        h = http.client.HTTPSConnection(host)
    else:
        h = http.client.HTTPConnection(host)

    h.putrequest('POST', selector)
    h.putheader('content-type', content_type)
    h.putheader('content-length', str(len(body)))
    h.endheaders()
    h.send(body)
    response = h.getresponse()
    return response.read()

def encode_multipart_formdata(fields, files):
    """
    fields is a sequence of (name, value) elements for regular form fields.
    files is a sequence of (name, filename, value) elements for
    data to be uploaded as files
    Return (content_type, body) ready for http.client connection instance
    """
    BOUNDARY_STR = '----------ThIs_Is_tHe_bouNdaRY_$'
    CRLF = bytes("\r\n","ASCII")
    L = []
    for (key, value) in fields:
        L.append(bytes("--" + BOUNDARY_STR,"ASCII"))
        L.append(bytes('Content-Disposition: form-data; name="%s"' % key,"ASCII"))
        L.append(b'')
        L.append(bytes(value,"ASCII"))
    for (key, filename, value) in files:
        L.append(bytes('--' + BOUNDARY_STR,"ASCII"))
        L.append(bytes('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename),"ASCII"))
        L.append(bytes('Content-Type: %s' % get_content_type(filename),"ASCII"))
        L.append(b'')
        L.append(value)
    L.append(bytes('--' + BOUNDARY_STR + '--',"ASCII"))
    L.append(b'')
    body = CRLF.join(L)
    content_type = 'multipart/form-data; boundary=' + BOUNDARY_STR
    return content_type, body

def get_content_type(filename):
    return mimetypes.guess_type(filename)[0] or 'application/octet-stream'