## SIR model using Javascript

[Link to Observable notebook](http://beta.observablehq.com/@epichef/deterministic-sir-model)

<script type="application/javascript">

function resizeIFrameToFitContent( iFrame ) {

    iFrame.width  = iFrame.contentWindow.document.body.scrollWidth;
    iFrame.height = iFrame.contentWindow.document.body.scrollHeight;
}

window.addEventListener('DOMContentLoaded', function(e) {

    var iFrame = document.getElementById( 'iFrame1' );
    resizeIFrameToFitContent( iFrame );

    // or, to resize all iframes:
    var iframes = document.querySelectorAll("iframe");
    for( var i = 0; i < iframes.length; i++) {
        resizeIFrameToFitContent( iframes[i] );
    }
} );

function importParentStyles() {
    var parentStyleSheets = parent.document.styleSheets;
    var cssString = "";
    for (var i = 0, count = parentStyleSheets.length; i < count; ++i) {
        if (parentStyleSheets[i].cssRules) {
            var cssRules = parentStyleSheets[i].cssRules;
            for (var j = 0, countJ = cssRules.length; j < countJ; ++j)
                cssString += cssRules[j].cssText;
        }
        else
            cssString += parentStyleSheets[i].cssText;  // IE8 and earlier
    }
    var style = document.createElement("style");
    style.type = "text/css";
    try {
        style.innerHTML = cssString;
    }
    catch (ex) {
        style.styleSheet.cssText = cssString;  // IE8 and earlier
    }
    document.getElementsByTagName("head")[0].appendChild(style);
}

</script>

<iframe src="../../observables/deterministic-sir-model/index.html" onload="this.width=screen.width;this.height=screen.height;" frameBorder="0" scrolling="no"></iframe>
