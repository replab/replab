var initCollapsedAdmonitions = () => {
  var admonitions = document.querySelectorAll('div.admonition');

  // Add the button to each admonition and hook up a callback to toggle visibility
  admonitions.forEach((item, index) => {
    var itemTitle = item.querySelector('.admonition-title');
    var admonitionID = `admonition-${index}`;
    item.setAttribute('id', admonitionID);
    var collapseButton = `<button id="button-${admonitionID}" class="admonition-button" data-target="${admonitionID}">+</button>`;
    itemTitle.insertAdjacentHTML('afterend', collapseButton);
    thisButton = document.getElementById(`button-${admonitionID}`);
    thisButton.addEventListener('click', toggleCollapsed);
  })
};

// This should simply add / remove the collapsed class and change the button text
var toggleCollapsed = (click) => {
  target = click.target.dataset['target']
  var admonition = document.getElementById(target);
  var button = admonition.querySelector('button.admonition-button');
  if (admonition.classList.contains("collapsed")) {
    admonition.classList.remove("collapsed");
    button.textContent = "-";
  } else {
    admonition.classList.add("collapsed");
    button.textContent = "+";
  }
}

// Helper function to run when the DOM is finished
const sphinxAdmonitionRunWhenDOMLoaded = cb => {
  if (document.readyState != 'loading') {
    cb()
  } else if (document.addEventListener) {
    document.addEventListener('DOMContentLoaded', cb)
  } else {
    document.attachEvent('onreadystatechange', function() {
      if (document.readyState == 'complete') cb()
    })
  }
}
sphinxAdmonitionRunWhenDOMLoaded(initCollapsedAdmonitions)
