var removeActive = function(){
	$("ul.navbar-nav").children().each(function() {
		$(this).removeClass('active');
	});	
};


var Router = Backbone.Router.extend({
	
	el: "#page-content", // selector to load page content

	routes: {
		'': 'home',
		'similarity(/)(:searchStr)': 'similarity',
		'drug/:id': 'drug',
		'clustergram': 'clustergram',
		'clusters': 'clusters',
		'downloads': 'downloads',
		'contribute': 'contribute'
	},

	home: function(){
		$(this.el).load("home.html", function() {
			removeActive();
			$("#browse").addClass("active");
		});
	},

	similarity: function(searchStr){
		$(this.el).load("similarity.html", function() {
			removeActive();
			$("#similarity").addClass('active');
			if (searchStr){
				doStrSearch(searchStr);	
			}
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
			$('#browse').addClass('active');	
		})
	},

	downloads: function() {
		$(this.el).load("downloads.html", function() {
			removeActive();
			$("#downloads").addClass('active');
		})
	},

	contribute: function() {
		$(this.el).load("contribute.html", function() {
			removeActive();
			$("#contribute").addClass('active');
		})
	},

});

var appRouter = new Router();
Backbone.history.start(); 

