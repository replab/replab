$(document).ready(function() {
    if (window.location.hash.length != 0) {
        base = window.location.hash.replace(/\./g, '\\.');
        base = base.replace(/\+/g, '\\+');
        admonitions = $('.admonition:has(' + base + ')');
	for (i = 0; i < admonitions.length; i ++) {
	    admonition = admonitions[i];
	    if (admonition.classList.contains("collapsed")) {
		admonition.classList.remove("collapsed");
	    }
	}
    }
    $('a').click(function(){
	target = $( this ).attr('href');
	pos = target.indexOf('#');
	if (pos != -1) {
	    hash = target.substr(pos);
            base = hash.replace(/\./g, '\\.');
            base = base.replace(/\+/g, '\\+');
            admonitions = $('.admonition:has(' + base + ')');
	    for (i = 0; i < admonitions.length; i ++) {
		admonition = admonitions[i];
		if (admonition.classList.contains("collapsed")) {
		    admonition.classList.remove("collapsed");
		}
	    }
	}
    });
});
