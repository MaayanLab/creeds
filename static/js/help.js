$(document).ready(function(){
	hljs.initHighlightingOnLoad();
	renderAllApiDocs();
})

function loadScript(scriptPath, obj){
	$.ajax({
		url: scriptPath,
		dataType: 'text',
		success: function(script){
			var suffix = scriptPath.split('.')[1];
			switch (suffix){
				case 'py':
					var scriptType = 'python';
					break;
				case 'json':
					var scriptType = 'json';
					break;
				default:
					var scriptType = 'r';
			}
			var codeBlock = $('<code>');
			codeBlock = codeBlock.html(script).addClass(scriptType);
			codeBlock = $('<pre>').append(codeBlock)
			obj.append(codeBlock);
		}
	});
}

function renderParams(data){
	var table = $('<table>').addClass('table parameters').append('<tbody>');
	var params = data["Parameters"]
	for (var i = 0; i < params.length; i++) {
		var param = params[i];
		var tr = $('<tr>');
		for (var j = 0; j < param.length; j++) {
			tr.append( $('<td>').html(param[j]) )
		};
		table.append(tr)
	};
	return table;
}

function renderScripts(data){
	var navTabs = $('<ul>').addClass('nav nav-tabs');
	var tabContents = $('<div>').addClass('tab-content');
	scriptPaths = data['Example code'];

	i = 0;
	for (var key in scriptPaths){
		var scriptPath = scriptPaths[key];
		var id = scriptPath.replace('.', '-');

		var li = $('<li>').attr('data-toggle', 'tab')
			.attr('role', 'presentation')
			.attr('data-target', '#'+id);
		li.append($('<a>').html(key));

		var tabPane = $('<div>').attr('id', id)
			.addClass('tab-pane fade');

		if (i === 0){
			li = li.addClass('active');
			tabPane = tabPane.addClass('active in');
		}

		navTabs = navTabs.append(li);
		tabContents = tabContents.append(tabPane);
		i++;
	}
	return {navTabs: navTabs, tabContents: tabContents};
}

function renderApiDoc(data){
	FIELDS = [
		'Method',
		'URL',
		'Returns',
		'Parameters',
		'Example code',
		'Example result',
		];
	var li = $('<li>').attr('id', data['id']);
	var h4 = $('<h4>').html(data['name']);
	var dl = $('<dl>').addClass('dl-horizontal');
	li = li.append(h4)
	
	for (var i = 0; i < FIELDS.length; i++) {
		var field = FIELDS[i];
		var dt = $('<dt>').html(field);
		var dd = $('<dd>')
		switch (field){
			case 'Parameters':
				dd.append(renderParams(data));
				break;
			case 'Example code':
				scriptDOMs = renderScripts(data)
				dd.append(scriptDOMs['navTabs']);
				dd.append(scriptDOMs['tabContents']);
				break;
			case 'Example result':
				dd.attr('id', data['Example result'].replace('.', '-'));
				break;
			default:
				dd.html(data[field]);
		}

		dl.append(dt);
		dl.append(dd);
	};

	li = li.append(dl);
	return li;
}

function renderAllApiDocs(){
	$.getJSON('api_docs.json', function(apiDocs){
		for (var i = 0; i < apiDocs.length; i++) {
			var li = renderApiDoc(apiDocs[i]);
			$('#api-docs').append(li);
			// append scripts after json got
			for (var key in apiDocs[i]['Example code']){
				var scriptPath = apiDocs[i]['Example code'][key];
				var id = scriptPath.replace('.', '-');
				loadScript(scriptPath, $('#' + id));
			}
			// append result
			var resultPath = apiDocs[i]['Example result'];
			var id = resultPath.replace('.', '-')
			loadScript(resultPath, $('#' + id))
			
		};
	});	
}