CREATE TABLE IF NOT EXISTS `hgnc` (
  `hgnc_id` varchar(16) NOT NULL,
  `approved_symbol` varchar(32) NOT NULL,
  `approved_name` text NOT NULL,
  `entrez_id` int(11) NOT NULL,
  PRIMARY KEY (`approved_symbol`)
);

CREATE TABLE IF NOT EXISTS `hgnc_synonyms` (
  `symbol` varchar(32) NOT NULL,
  `synonym` text NOT NULL
);

CREATE TABLE IF NOT EXISTS `homologene` (
  `homologene_id` int(11) NOT NULL,
  `hgnc_symbol` varchar(32) NOT NULL,
  `mgi_symbol` varchar(32) NOT NULL,
  UNIQUE (`homologene_id`)
);

CREATE TABLE IF NOT EXISTS `mgi` (
  `mgi_id` varchar(32) NOT NULL,
  `mgi_symbol` varchar(32) NOT NULL,
  `mgi_name` text,
  `type` varchar(64) NOT NULL,
  `entrez_id` int(11) NOT NULL,
  PRIMARY KEY (`mgi_id`)
);

CREATE TABLE IF NOT EXISTS `mgi_synonyms` (
  `mgi_symbol` varchar(32) NOT NULL,
  `mgi_synonym` text NOT NULL
);
