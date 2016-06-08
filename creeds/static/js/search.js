ENTER_POINT = '/CREEDS';
var exampleGenes = {
	upGenes: ['KIAA0907','KDM5A','CDC25A','EGR1','GADD45B','RELB','TERF2IP','SMNDC1','TICAM1','NFKB2','RGS2','NCOA3','ICAM1','TEX10','CNOT4','ARID4B','CLPX','CHIC2','CXCL2','FBXO11','MTF2','CDK2','DNTTIP2','GADD45A','GOLT1B','POLR2K','NFKBIE','GABPB1','ECD','PHKG2','RAD9A','NET1','KIAA0753','EZH2','NRAS','ATP6V0B','CDK7','CCNH','SENP6','TIPARP','FOS','ARPP19','TFAP2A','KDM5B','NPC1','TP53BP2','NUSAP1'],
	dnGenes: ['SCCPDH','KIF20A','FZD7','USP22','PIP4K2B','CRYZ','GNB5','EIF4EBP1','PHGDH','RRAGA','SLC25A46','RPA1','HADH','DAG1','RPIA','P4HA2','MACF1','TMEM97','MPZL1','PSMG1','PLK1','SLC37A4','GLRX','CBR3','PRSS23','NUDCD3','CDC20','KIAA0528','NIPSNAP1','TRAM2','STUB1','DERA','MTHFD2','BLVRA','IARS2','LIPA','PGM1','CNDP2','BNIP3','CTSL1','CDC25B','HSPA8','EPRS','PAX8','SACM1L','HOXA5','TLE1','PYGL','TUBB6','LOXL1']
};


$(document).ready(function(){
	// autocomplete for searching signatures
	$.getJSON('data/autoCompleteList.json', function(signatureInfo){
		// console.log('autoCompleteList got')
		// console.log(signatureInfo)

		var autoCompleteGenes = new Bloodhound({
			datumTokenizer: Bloodhound.tokenizers.whitespace,
			queryTokenizer: Bloodhound.tokenizers.whitespace,
			local: signatureInfo.gene
		});

		var autoCompleteDzs = new Bloodhound({
			datumTokenizer: Bloodhound.tokenizers.whitespace,
			queryTokenizer: Bloodhound.tokenizers.whitespace,
			local: signatureInfo.dz
		});

		var autoCompleteDrugs = new Bloodhound({
			datumTokenizer: Bloodhound.tokenizers.whitespace,
			queryTokenizer: Bloodhound.tokenizers.whitespace,
			local: signatureInfo.drug
		});

		$("#sigNameInput").typeahead({
			hint: true,
			highlight: true,
			minLength: 2
		},
		{
			name: 'genes',
			source: autoCompleteGenes,
			templates: {
				header: '<h4 class="catagory-header">Genes</h4>'
			}
		},
		{
			name: 'dzs',
			source: autoCompleteDzs,
			templates: {
				header: '<h4 class="catagory-header">Diseases</h4>'
			}
		},
		{
			name: 'drugs',
			source: autoCompleteDrugs,
			templates: {
				header: '<h4 class="catagory-header">Drugs</h4>'
			}
		}
		);
	});

	
	$("#sigNameBtn").click(function(){
		var searchStr = $("#sigNameInput").val();
		var baseUrl = window.location.origin + window.location.pathname + '#similarity' ;
		var redirectUrl = baseUrl + '/' + searchStr
		window.location.assign(redirectUrl);
	});

	$(".exampleTerm").click(function(){
		var term = $(this).text();
		$("#sigNameInput").val(term);
	})

	$("#geneSearchEgBtn").click(function(){
		$('#upGenes').val(exampleGenes.upGenes.join('\n'));
		$('#dnGenes').val(exampleGenes.dnGenes.join('\n'));
	});

	// search sigs using up/dn genes
	$(".geneSearchBtn").click(function(e){
		e.preventDefault();
		var upGenes = $('#upGenes').val();
		var dnGenes = $('#dnGenes').val();
		var dbVersion = [];
		$.each($(".dbVersion input:checked"), function(){
			dbVersion.push($(this).val());
		});

		if (dbVersion.length === 0) {
			alert('Please select at least one Database version');
		} 

		else if (upGenes != '' && dnGenes !='') {
	        var direction = $(this).attr('direction');
			var postPayLoad = {"up_genes": upGenes.split('\n'), "dn_genes": dnGenes.split('\n'), 
				"direction":direction, "db_version": dbVersion, "client": "webapp"};
			$.ajax({
				type: 'POST',
				url: ENTER_POINT+'/search',
				contentType : 'application/json',
				data: JSON.stringify(postPayLoad),
				dataType: 'json',
				success: function(id){
					// redirect to result page
					var baseUrl = window.location.origin + window.location.pathname + '#query';
					var redirectUrl = baseUrl + '/' + id;
					window.location.assign(redirectUrl);
				},				

			});

		} else {
			alert('Please fill genes');
		};

	})
});


function submitSelectedBtn (){
	$("#submitSelectedBtn").remove();
	var btn = $("<a>");
	btn.text("Visualize results as a clustergram");
	btn.addClass('btn btn-info');
	btn.attr('href', '#clustergram');
	btn.attr('id', 'submitSelectedBtn');
	btn.click(function(){ // click to post genes and uids to clustergram vis
		
		var upGenes = $('#upGenes').val().split('\n');
		var dnGenes = $('#dnGenes').val().split('\n');		
		var genes = upGenes.concat(dnGenes);

		postPayLoad = {"genes": genes, "ids": uids}; // global

	});

	$("#searchResults").append(btn);

}

