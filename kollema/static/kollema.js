function closeModal() {
    for (i = 0; i < arguments.length; i++) {
        var modal = document.getElementById(arguments[i])
        modal.style.display='none'
    }
}

function openModal() {
    for (var i = 0; i < arguments.length; i++) {
        var modal = document.getElementById(arguments[i]);
        modal.style.display='block';
        ch = modal.getElementsByTagName('input');
        for (var j=0; j < ch.length; j++) {
            if ( ch[j].getAttribute("type") == "text") {
                ch[j].value='';
            }
        }
    }
}
