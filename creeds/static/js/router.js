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
		'query(/)(:id)': 'query',
		'clustergram': 'clustergram',
		'bubble': 'bubble',
		'clusters': 'clusters',
		'downloads': 'downloads',
		'help(/)(:target)': 'help'
	},

	home: function(){
		$(this.el).load("landing.html", function() {
			removeActive();
		});
	},

	bubble: function(){
		$(this.el).load("home.html", function() {
			removeActive();
			$("#browse").addClass("active");
		});		
	},

	similarity: function(searchStr){
		// display string search result
		$(this.el).load("similarity.html", function() {
			removeActive();
			$("#similarity").addClass('active');
			if (searchStr){
				doStrSearch(searchStr);	
			}
		});
	},

	query: function(id){
		// display genes query result
		$(this.el).load("similarity.html", function() {
			removeActive();
			$("#similarity").addClass("active");
			if (id){
				doQuery(id);
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

	help: function(target) {
		if ($("#help-page").length){ // current page is help page
			if (target){
				scrollTo(target);
			}
		}else {
			$(this.el).load("help.html", function() {
				removeActive();
				$("#help").addClass('active');
				$('html').on('custom', function(e, eventData){
					console.log(eventData)
					if (eventData === "api3_result.json"){
						// the last code block loaded
						if (target){
							scrollTo(target);
						}
					}
				});
			});			
		}
	},

});

var appRouter = new Router();
Backbone.history.start(); 

