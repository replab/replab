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
});
