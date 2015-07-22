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
		'clustergram': 'clustergram',
		'clusters': 'clusters'
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
			removeActive();
			$("#browse").addClass('active');			
		});
	},

	clustergram: function(){
		$(this.el).load("clustergram.html", function() {
			removeActive();
		})
	},
	
	clusters: function() {
		$(this.el).load("clusters.html", function() {
			removeActive();
			$('#clusters').addClass('active');	
		})
	},

});

var appRouter = new Router();
Backbone.history.start(); 