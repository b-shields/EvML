$(document).ready(function(){
  $.get("http://ipinfo.io", function(response) {
    Shiny.onInputChange("getIP", response);
  }, "json");
});