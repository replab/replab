var processAnchors = () => {
    if (window.location.hash.length != 0) {
        base = window.location.hash.replace(/\./g, '\\.');
        base = base.replace(/\+/g, '\\+');
	console.log(base);
        admonitions = $('.admonition:has(' + base + ')');
	console.log(admonitions);
	for (i = 0; i < admonitions.length; i ++) {
	    admonition = admonitions[i];
	    console.log(admonition.classList);
	    if (admonition.classList.contains("toggle-hidden")) {
		admonition.classList.remove("toggle-hidden");
	    }
	    buttons = $(admonition).children(".toggle-button");
	    console.log(buttons);
	    if (buttons.length == 1 && buttons[0].classList.contains("toggle-button-hidden"))
		buttons[0].classList.remove("toggle-button-hidden")
	}
    }
};

var watchAnchorClicks = () => {
    $(document).ready(function() {
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
		    if (admonition.classList.contains("toggle-hidden")) {
			admonition.classList.remove("toggle-hidden");
		    }
		}
		buttons = $(admonition).children(".toggle-button");
		if (buttons.length == 1 && buttons[0].classList.contains("toggle-button-hidden"))
		    buttons[0].classList.remove("toggle-button-hidden")
	    }
	});
    });
};

sphinxToggleRunWhenDOMLoaded(processAnchors);
sphinxToggleRunWhenDOMLoaded(watchAnchorClicks);
