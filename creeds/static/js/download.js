var FileMeta = Backbone.Model.extend({
	defaults: {
		category: '',
		num_sigs: 0,
		filenames: {},
	},
});

// View for FileMeta Model
var FileMetaView = Backbone.View.extend({
	model: FileMeta,
	tagName: 'tr',

	template: _.template($('#fileTemplate').html()),

	initialize: function(){
		this.render();
	},
	render: function(){
		this.$el.html( this.template(this.model.toJSON()));
		return this;
	},

});

var FileMetas = Backbone.Collection.extend({
	model: FileMeta,
	url: ENTER_POINT + '/download',
});

// View for Collection of FileMetas
var FileMetasView = Backbone.View.extend({
	tagName: 'tbody',
	initialize: function(){
		this.listenTo(this.collection, 'sync', this.render)

		this.collection.fetch();
	},

	render: function(){
		this.collection.each(function(fileMeta){
			var fileMetaView = new FileMetaView({model: fileMeta});
			this.$el.append(fileMetaView.el);
		}, this);
		return this;
	}
});


// Init instances
var fileMetas = new FileMetas();
var fileMetasView = new FileMetasView({ collection: fileMetas });
// Render view of the collection
$('table').append(fileMetasView.render().el);
