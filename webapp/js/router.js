var removeActive = function(){
	$("ul.navbar-nav").children().each(function() {
		$(this).removeClass('active');
	});	
};


var Router = Backbone.Router.extend({
	
	el: "#page-content", // selector to load page content

	routes: {
		'': 'home',
		'similarity': 'similarity',
		'drug/:id': 'drug',
	},

	home: function(){
		$(this.el).load("home.html", function() {
			removeActive();
			$("#home").addClass("active");
		});
	},

	similarity: function(){
		$(this.el).load("similarity.html", function() {
			removeActive();
			$("#similarity").addClass('active');
		});
	},

	drug: function(id){
		$(this.el).load("drug_profile.html", function() {
			hideTour();
			removeActive();
			$("#browse").addClass('active');			
		});
	},

});

var appRouter = new Router();
Backbone.history.start(); 