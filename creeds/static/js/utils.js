function doStrSearch(searchStr){
	// wrapper for send GET request to string search API and displayStrSearchResult
	$.getJSON(ENTER_POINT+'/search', {q: searchStr}, function(results){
		displayStrSearchResult(results);
	});

}

function displayStrSearchResult(results){
	d3.select("#searchResultsInner").remove()
	d3.select("#searchResults").append('div')
		.attr('id', 'searchResultsInner')

	d3.select("#searchResultsInner").append('div')
		.attr('class','alert alert-success')
		.text(results.length + ' signatures found');

	for (var i = results.length - 1; i >= 0; i--) {
		var info = results[i];

		domId = info['id'].replace(':','-')
		d3.select("#searchResultsInner")
			.append("div")
			.attr('id', domId)
			.attr('class', 'well')
			// .style("height", '600px')
			// .style("overflow", "auto")
		
		var div = d3.select("#" + domId); // the container to put node info
		
		// display the meta data of the signature in the dl
		var dl = div.append("dl")
			.attr('class', 'dl-horizontal')

		var prefix = info['id'].split(':')[0]
		switch (prefix) {
			case 'gene':
				var specificKeys = ['hs_gene_symbol', 'mm_gene_symbol', 'pert_type'];
				var label = 'single gene perturbation';
				var color = '#9467bd';
				break;
			case 'dz':
				var specificKeys = ['disease_name', 'umls_cui', {'do_id':'http://disease-ontology.org/term/'}];
				var label = 'disease';
				var color = '#1f77b4';
				break;
			case 'drug':
				var specificKeys = ['drug_name', {'drugbank_id':'http://www.drugbank.ca/drugs/'}, {'pubchem_cid':'https://pubchem.ncbi.nlm.nih.gov/compound/'}, 'smiles'];
				var label = 'single drug perturbation';
				var color = '#ff7f0e';
				break;
		}


		dl.append('dt').text('type');
		dl.append('dd').append('span')
			.attr('class', 'label')
			.style('background-color', color)
			.text(label);
		
		var profileUrl = ENTER_POINT + '/api?id=' + info['id'];
		dl.append('dt').text('id');
		dl.append('dd').append('a')
			.attr('href', profileUrl)
			.attr('target', '_blank')
			.text(info['id']);

		var keys = [{'geo_id': 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='}, 
			'ctrl_ids', 'pert_ids', {'platform': 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='}, 
			'organism', 'cell_type'] // the generic keys in info object
		var allKeys = specificKeys.concat(keys)
	    for (var ii = 0; ii < allKeys.length; ii++) {
	    	var key = allKeys[ii]
	    	if(typeof(key) === 'string'){
		    	dl.append('dt').text(key);
		    	var value = info[key];
		    	if (typeof(value) !== 'string' && value !== undefined){
		    		if (_.isObject(value[0])){ // array of objects for v2.0 signaturess
		    			var texts = _.map(info[key], function(d){ return d['name']; });
		    		} else {
		    			var texts = value;
		    		}
					dl.append('dd').text(texts.join(', '));
		    	} else { 
		    		dl.append('dd').text(info[key]);
		    	}
		    	
	    	} else { // add hyperlinks for key that is object
	    		var field = Object.keys(key)[0];
	    		var url = key[field] + info[field];
	    		dl.append('dt').text(field);
	    		dl.append('dd').append('a').attr('href', url)
	    			.attr('target', '_blank')
	    			.text(info[field]);
	    	}
	    }
	};
	scrollToResult();    	
}

function scrollToResult() {
	$('html, body').animate({
		scrollTop: $("#searchResults").offset().top
	}, 500);
}

function scrollTo(selector) {
	$('html, body').animate({
		scrollTop: $(selector).offset().top
	}, 500);
}

function doQuery(id){
	// wrapper for send GET request to /result endpoint and displayGeneSearchResult
	$.blockUI({ css: { 
		            border: 'none', 
		            padding: '15px', 
		            backgroundColor: '#000', 
		            '-webkit-border-radius': '10px', 
		            '-moz-border-radius': '10px', 
		            opacity: .5, 
		            color: '#fff' 
		        } });

	$.getJSON(ENTER_POINT+'/result', {id: id}, function(results){
		displayGeneSearchResult(results);
		$.unblockUI();
		$('[data-toggle="tooltip"]').tooltip();
	});
}


function displayGeneSearchResult (results) { // to display results using gene list search
	d3.select('#searchResultsInner').remove();
	d3.select("#searchResults").append('div')
		.attr('id', 'searchResultsInner');	

	headerLabels = ['ID', 'Name', 'GEO ID', 'Signed Jaccard Index ']
	var table = $('<table>').addClass('table table-striped table-hover');

	var thead = $('<thead>');
	var tr = $('<tr>')

	$.each(headerLabels,function(i, headerLabel) {
		tr.append('<th>' + headerLabel + '</th>');
	});

	thead.append(tr)
	table.append(thead);

	$("#searchResultsInner").append(table);

	uids = [];

	var columnsDef = [
		{ 
			"data": "id",
			"render": function(data, type, full, meta){
				var profileUrl = ENTER_POINT + '/api?id=' + data;
				if (data.startsWith('dz:')){
					// prettify the results for displaying
					data = 'disease:' + data.split(':')[1];
				}
				return '<a target="_blank" href="'+profileUrl+'">'+data+'</a>';
			}
		},
	    {
	    	"data": "name", 
			"render": function(data, type, full, meta){
				if (_.isArray(data[0])){
					html = '';
					for (var i = 0; i < data[0].length; i++) {
						var url = data[1][i];
						var name = data[0][i];
						var a = '<a target="_blank" href="'+url+'">'+name+'</a><br>';
						html += a;
					};
					html += ''
					return html;

				}else {
					return '<a target="_blank" href="'+data[1]+'">'+data[0]+'</a>';	
				}
				
			}
		},
	    { 
	    	"data": "geo_id",
	    	"render": function(data, type, full, meta){
	    		return '<a target="_blank" href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+data+'">'+data+'</a>';
	    	}
	    },
		{ "data": "signed_jaccard" },
	];

	rctable = table.dataTable({
		"data": results,
		"columns": columnsDef,
		"bSort": false,
		"deferRender": true,
		"sDom": 'T<"clear">lfrtip',
		"oTableTools": {
        	// to save table
        	"sSwfPath": "js/swf/copy_csv_xls_pdf.swf",
        	"aButtons": [
        	"copy",
        	"print",
        	{
        		"sExtends":    "collection",
        		"sButtonText": 'Save <span class="caret" />',
        		"aButtons":    [ "csv", "xls", "pdf" ]
        	}
        	],
            // to make rows selectable
            // "sRowSelect": "multi",
            // "fnRowSelected": function ( row ) {
            // 	var rData = rctable.fnGetData(row);
            // 	uids.push(rData.id);
            // },
            // "fnRowDeselected": function (row) {
            // 	var rData = rctable.fnGetData(row);
            // 	// remove from uids
            // 	var index = uids.indexOf(rData.id);
            // 	if (index > -1) {
            // 		uids.splice(index, 1);
            // 	}
            // },
        },
        "fnInitComplete": function(oSettings, json){ //DataTables has finished its initialisation.
        	var infoIcon = $('<i title="Leverages the effect of directions between a pair of gene expression signatures. It has a range of [-1,1] where 1 represent completely same signatures and -1 represent signatures of reverse effects and 0 represent unrelated signatures." data-toggle="tooltip" data-placement="right" class="glyphicon glyphicon-info-sign"></i>');
        	$('th').last().append(infoIcon);
        },
    });
	scrollToResult();
}
