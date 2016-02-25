/************************************************************************
@file GetConceptCovariates.sql

Copyright 2015 Observational Health Data Sciences and Informatics

This file is part of CohortMethod

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
************************************************************************/

{DEFAULT @cdm_database = 'CDM4_SIM' } 
{DEFAULT @cdm_version == '4'}
{DEFAULT @cohort_table = '#cohort_person'}
{DEFAULT @row_id_field = 'row_id'}
{DEFAULT @cohort_definition_id = 'cohort_concept_id'} 
{DEFAULT @concept_class_id = 'concept_class'} 
{DEFAULT @measurement = 'observation'} 
{DEFAULT @use_covariate_demographics = TRUE}
{DEFAULT @use_measurement = TRUE}

USE @cdm_database;

IF OBJECT_ID('tempdb..#cov', 'U') IS NOT NULL
	DROP TABLE #cov;

IF OBJECT_ID('tempdb..#cov_ref', 'U') IS NOT NULL
	DROP TABLE #cov_ref;

CREATE TABLE #cov_ref (
	covariate_id BIGINT,
	covariate_name VARCHAR(512),
	analysis_id INT,
	concept_id INT
	);
	


/**************************
***************************
DEMOGRAPHICS
***************************
**************************/
{@use_covariate_demographics} ? {

--gender
SELECT cp1.@row_id_field AS row_id,
	gender_concept_id AS covariate_id,
	1 AS covariate_value
INTO #cov_gender
FROM @cohort_table cp1
INNER JOIN person p1
	ON cp1.subject_id = p1.person_id
WHERE p1.gender_concept_id IN (
		SELECT concept_id
		FROM concept
		WHERE LOWER(@concept_class_id) = 'gender'
		);


INSERT INTO #cov_ref (
  covariate_id,
	covariate_name,
	analysis_id,
	concept_id
	)
SELECT p1.covariate_id,
	'Gender = ' +
    CASE WHEN v1.concept_name IS NOT NULL
			THEN v1.concept_name
		ELSE 'Unknown invalid concept'
		END AS covariate_name,
	2 AS analysis_id,
	p1.covariate_id AS concept_id
FROM (SELECT distinct covariate_id FROM #cov_gender) p1
LEFT JOIN (
	SELECT concept_id,
		concept_name
	FROM concept
	WHERE LOWER(@concept_class_id) = 'gender'
	) v1
	ON p1.covariate_id = v1.concept_id;

--race
SELECT cp1.@row_id_field AS row_id,
	race_concept_id AS covariate_id,
	1 AS covariate_value
  INTO #cov_race
FROM @cohort_table cp1
INNER JOIN person p1
	ON cp1.subject_id = p1.person_id
WHERE p1.race_concept_id IN (
		SELECT concept_id
		FROM concept
		WHERE LOWER(@concept_class_id) = 'race'
		);


INSERT INTO #cov_ref (
  covariate_id,
  covariate_name,
	analysis_id,
	concept_id
	)
SELECT p1.covariate_id,
	'Race = ' + CASE WHEN v1.concept_name IS NOT NULL
  		THEN v1.concept_name
		ELSE 'Unknown invalid concept'
		END  AS covariate_name,
	3 AS analysis_id,
	p1.covariate_id AS concept_id
FROM (SELECT distinct covariate_id FROM #cov_race) p1
LEFT JOIN (
	SELECT concept_id,
		concept_name
	FROM concept
	WHERE LOWER(@concept_class_id) = 'race'
	) v1
	ON p1.covariate_id = v1.concept_id;

--ethnicity
SELECT cp1.@row_id_field AS row_id,
	ethnicity_concept_id AS covariate_id,
	1 AS covariate_value
  INTO #cov_ethnicity
FROM @cohort_table cp1
INNER JOIN person p1
	ON cp1.subject_id = p1.person_id
WHERE p1.ethnicity_concept_id IN (
		SELECT concept_id
		FROM concept
		WHERE LOWER(@concept_class_id) = 'ethnicity'
		);



INSERT INTO #cov_ref (
  covariate_id,
  covariate_name,
  analysis_id,
	concept_id
	)
SELECT p1.covariate_id,
	'Ethnicity = ' + CASE WHEN v1.concept_name IS NOT NULL
  		THEN v1.concept_name
		ELSE 'Unknown invalid concept'
		END  AS covariate_name,
	4 AS analysis_id,
	p1.covariate_id AS concept_id
FROM (SELECT distinct covariate_id FROM #cov_ethnicity) p1
LEFT JOIN (
	SELECT concept_id,
		concept_name
	FROM concept
	WHERE LOWER(@concept_class_id) = 'ethnicity'
	) v1
	ON p1.covariate_id = v1.concept_id;


-- age day:
SELECT cp1.@row_id_field AS row_id,
  0 AS covariate_id,
	datediff(day, datefromparts(p1.YEAR_OF_BIRTH,isnull(p1.MONTH_OF_BIRTH,1),1), cp1.cohort_start_date) AS covariate_value
    INTO #cov_age_day
FROM @cohort_table cp1
INNER JOIN person p1
	ON cp1.subject_id = p1.person_id
WHERE (YEAR(cp1.cohort_start_date) - p1.YEAR_OF_BIRTH) >= 0
	AND (YEAR(cp1.cohort_start_date) - p1.YEAR_OF_BIRTH) < 100;


INSERT INTO #cov_ref (
  covariate_id,
	covariate_name,
	analysis_id,
	concept_id
	)
SELECT p1.covariate_id,
	'Age in days'  AS covariate_name,
	4 AS analysis_id,
	0 AS concept_id
FROM (select distinct covariate_id FROM #cov_age_day) p1
;


--age group
SELECT cp1.@row_id_field AS row_id,
	FLOOR((YEAR(cp1.cohort_start_date) - p1.YEAR_OF_BIRTH) / 5) + 10 AS covariate_id,
	1 AS covariate_value
    INTO #cov_age
FROM @cohort_table cp1
INNER JOIN person p1
	ON cp1.subject_id = p1.person_id
WHERE (YEAR(cp1.cohort_start_date) - p1.YEAR_OF_BIRTH) >= 0
	AND (YEAR(cp1.cohort_start_date) - p1.YEAR_OF_BIRTH) < 100;




INSERT INTO #cov_ref (
  covariate_id,
	covariate_name,
	analysis_id,
	concept_id
	)
SELECT p1.covariate_id,
	'Age group: ' + CAST((covariate_id-10)*5 AS VARCHAR) + '-' + CAST((covariate_id-10+1)*5-1 AS VARCHAR)  AS covariate_name,
	4 AS analysis_id,
	0 AS concept_id
FROM (select distinct covariate_id FROM #cov_age) p1
;



--index month

SELECT cp1.@row_id_field AS row_id,
	MONTH(cohort_start_date) + 40 AS covariate_id,
	1 AS covariate_value
    INTO #cov_month
FROM @cohort_table cp1;


INSERT INTO #cov_ref (
  covariate_id,
  covariate_name,
  analysis_id,
	concept_id
	)
SELECT p1.covariate_id,
	'Index month: ' + CAST(covariate_id-40 AS VARCHAR)  AS covariate_name,
	6 AS analysis_id,
	0 AS concept_id
FROM (select distinct covariate_id FROM #cov_month) p1
;

}

/**************************
***************************
RISK FACTOR CONCEPTS
***************************
**************************/

IF OBJECT_ID('tempdb..#risk_factors_all', 'U') IS NOT NULL
	DROP TABLE #risk_factors_all;
	
	select * into #risk_factors_all from 
	(select a.covariate_id, b.descendant_concept_id concept_id, a.priorTime
	from #risk_factors a inner join concept_ancestor b
	on a.concept_id=b.ancestor_concept_id
	union select covariate_id, concept_id, priorTime from #risk_factors) temp
	;

-- 1) get the conditions
IF OBJECT_ID('tempdb..#rf_temp_con', 'U') IS NOT NULL
	DROP TABLE #rf_temp_con;
	
	select distinct a.@row_id_field, a.subject_id, c.covariate_id, b.condition_start_date start_date, 
	b.condition_concept_id concept_id 
	into #rf_temp_con
	from
	@cohort_table a inner join condition_occurrence b 
	on a.subject_id=b.person_id 
	inner join #risk_factors_all c
	on b.condition_concept_id=c.concept_id
	inner join observation_period d
	on a.subject_id=d.person_id 
	where datediff(day, b.condition_start_date, a.cohort_start_date) between 1 and c.priorTime
	and b.condition_start_date between d.observation_period_start_date and d.observation_period_end_date
	;
-- 2) get the drugs
IF OBJECT_ID('tempdb..#rf_temp_drug', 'U') IS NOT NULL
	DROP TABLE #rf_temp_drug;
	
	select distinct a.@row_id_field, a.subject_id, c.covariate_id, b.drug_exposure_start_date start_date, 
	b.drug_concept_id concept_id 
	into #rf_temp_drug
	from
	@cohort_table a inner join drug_exposure b 
	on a.subject_id=b.person_id 
	inner join #risk_factors_all c
	on b.drug_concept_id=c.concept_id
	inner join observation_period d
	on a.subject_id=d.person_id 
	where datediff(day, b.drug_exposure_start_date, a.cohort_start_date) between 1 and c.priorTime
	and b.drug_exposure_start_date between d.observation_period_start_date and d.observation_period_end_date
		;
-- 3) get the procedures
IF OBJECT_ID('tempdb..#rf_temp_pro', 'U') IS NOT NULL
	DROP TABLE #rf_temp_pro;
	
	select distinct a.@row_id_field, a.subject_id, c.covariate_id, b.procedure_date start_date, 
	b.procedure_concept_id concept_id 
	into #rf_temp_pro
	from
	@cohort_table a inner join procedure_occurrence b 
	on a.subject_id=b.person_id 
	inner join #risk_factors_all c
	on b.procedure_concept_id=c.concept_id
	inner join observation_period d
	on a.subject_id=d.person_id 
	where datediff(day, b.procedure_date, a.cohort_start_date) between 1 and c.priorTime
	and b.procedure_date between d.observation_period_start_date and d.observation_period_end_date
		;
-- 4) get the observation
IF OBJECT_ID('tempdb..#rf_temp_obs', 'U') IS NOT NULL
	DROP TABLE #rf_temp_obs;
	
	select distinct a.@row_id_field, a.subject_id, c.covariate_id, b.observation_date start_date, 
	b.observation_concept_id concept_id 
	into #rf_temp_obs
	from
	@cohort_table a inner join observation b 
	on a.subject_id=b.person_id 
	inner join #risk_factors_all c
	on b.observation_concept_id=c.concept_id
	inner join observation_period d
	on a.subject_id=d.person_id 
	where datediff(day, b.observation_date, a.cohort_start_date) between 1 and c.priorTime
	and b.observation_date between d.observation_period_start_date and d.observation_period_end_date;
			
{@use_measurement} ? {
IF OBJECT_ID('tempdb..#rf_temp_mea', 'U') IS NOT NULL
	DROP TABLE #rf_temp_mea;
	
	select distinct a.@row_id_field, a.subject_id, c.covariate_id, b.measurement_date start_date, 
	b.measurement_concept_id concept_id 
	into #rf_temp_mea
	from
	@cohort_table a inner join measurement b 
	on a.subject_id=b.person_id 
	inner join #risk_factors_all c
	on b.measurement_concept_id=c.concept_id
	inner join observation_period d
	on a.subject_id=d.person_id 
	where datediff(day, b.measurement_date, a.cohort_start_date) between 1 and c.priorTime
	and b.measurement_date between d.observation_period_start_date and d.observation_period_end_date;
		
}

--rf_all table
IF OBJECT_ID('tempdb..#rf_all', 'U') IS NOT NULL
	DROP TABLE #rf_all;
-- put all the tables together to get covariate values:
select row_id, covariate_id, count(distinct start_date) covariate_value
INTO #rf_all
from 
(select @row_id_field row_id, covariate_id, start_date
from #rf_temp_con 
union 
select @row_id_field row_id, covariate_id,  start_date
from #rf_temp_drug 
union 
select @row_id_field row_id, covariate_id,  start_date 
from #rf_temp_obs 
union 
select @row_id_field row_id, covariate_id,  start_date
from #rf_temp_pro
{@use_measurement} ? {
union 
select @row_id_field row_id, covariate_id, start_date 
from #rf_temp_mea}
 ) temp
 group by row_id, covariate_id;
	
-- add the info into the covariate ref
INSERT INTO #cov_ref (
  covariate_id,
  covariate_name,
  analysis_id,
	concept_id
	)
SELECT distinct covariate_id,
	covariate_name  AS covariate_name,
	-1 AS analysis_id,
	0 as concept_id
FROM #risk_factors;	
	
	
	
-- combine rf_all and demo if wanted into #cov:

select * into #cov
from
(select * from #rf_all
{@use_covariate_demographics} ? {
union select * from #cov_gender
union select * from #cov_race
union select * from #cov_ethnicity
union select * from #cov_age
union select * from #cov_age_day
union select * from #cov_month
}) temp;
	


IF OBJECT_ID('tempdb..#cov_exposure', 'U') IS NOT NULL
  DROP TABLE #cov_exposure;
IF OBJECT_ID('tempdb..#cov_gender', 'U') IS NOT NULL
  DROP TABLE #cov_gender;
IF OBJECT_ID('tempdb..#cov_race', 'U') IS NOT NULL
  DROP TABLE #cov_race;
IF OBJECT_ID('tempdb..#cov_ethnicity', 'U') IS NOT NULL
  DROP TABLE #cov_ethnicity;
IF OBJECT_ID('tempdb..#cov_age', 'U') IS NOT NULL
  DROP TABLE #cov_age;
IF OBJECT_ID('tempdb..#cov_age_day', 'U') IS NOT NULL
  DROP TABLE #cov_age_day;
IF OBJECT_ID('tempdb..#cov_month', 'U') IS NOT NULL
  DROP TABLE #cov_month;
IF OBJECT_ID('tempdb..#rf_temp_con', 'U') IS NOT NULL
  DROP TABLE #rf_temp_con;
IF OBJECT_ID('tempdb..#rf_temp_drug', 'U') IS NOT NULL
  DROP TABLE #rf_temp_drug;
IF OBJECT_ID('tempdb..#rf_temp_obs', 'U') IS NOT NULL
  DROP TABLE #rf_temp_obs;
IF OBJECT_ID('tempdb..#rf_temp_pro', 'U') IS NOT NULL
  DROP TABLE #rf_temp_pro;
{@use_measurement}?{
IF OBJECT_ID('tempdb..#rf_temp_mea', 'U') IS NOT NULL
  DROP TABLE #rf_temp_mea;
}
IF OBJECT_ID('tempdb..#rf_all', 'U') IS NOT NULL
  DROP TABLE #rf_all;


