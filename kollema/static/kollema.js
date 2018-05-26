function closeModal() {
    for (i = 0; i < arguments.length; i++) {
        //var modal = document.getElementById(arguments[i])
        //modal.style.display='none'
        var target = '#' + arguments[i];
        $(target).hide();
    }
}

function openModal(source, element_id) {
    //alert('n ' + arguments.length+' '+arguments[0])
    var url = 'static/' + source + '.html';
    var target = '#' + element_id;
    //alert('source '+source+'     target '+target+arguments[0])
    $(target).load( url );
    //var modal = document.getElementById(element_id);
    //modal.style.display='block';
    $(target).show();
    console.log('showing '+ target)

    for (var i = 2; i < arguments.length; i++) {
        var t = '#' + arguments[i];
        //alert('id ' + target);
        $(t).show();
        console.log('showing '+ t)
        $(t + ' input').value='';
//        ch = modal.getElementsByTagName('input');
//        for (var j=0; j < ch.length; j++) {
//            // reset text in all input boxes
//            if ( ch[j].getAttribute("type") == "text") {
//                ch[j].value='';
//            }
//        }
    }
}

function openModalUrl(url, target) {
/*--------------------------------------------------------------------------------------------------
Load a url into a target, typically a modal dialog
url: address of target html document, possibly containing id selector

example:
<a class="menubutton" onclick="openModalUrl('static/transcript.html #load', '#reg')">

:param: url, address of html page, possibly containing sub-selector
:param: target, html id of page element
--------------------------------------------------------------------------------------------------*/
    $(target).load( url );
    $(target+' input').value='';
    $(target).show();
    //console.log('showing '+ target);
}

function getSessionVar( key ) {
    var request = $.ajax({'url': '/getSessionVar?key=' + key});
    alert('getSessionVar '+request[key]);
    return request[key];
}
