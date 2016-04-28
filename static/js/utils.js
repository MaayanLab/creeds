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
		    	dl.append('dd').text(info[key]);	    		
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
