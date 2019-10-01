package pt.uminho.ceb.biosystems.merlin.interpro.services;

import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;
import java.util.Observable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.Xref;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.interpro.InterProResult;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.interpro.InterProResultsList;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.interpro.Location;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.interpro.Model;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.InterproStatus;
import pt.uminho.ceb.biosystems.merlin.services.interpro.InterproServices;

/**
 * @author Oscar Dias
 *
 */
public class InterProDataLoadingServices extends Observable implements Runnable {

	private Statement statement;
	private AtomicBoolean cancel;
	private Map<String, InterProResultsList> interProResultsLists;
	private AtomicInteger datum;
	private ConcurrentLinkedQueue<String> list;
	private ConcurrentHashMap<String, Integer> auxiliaryMap;
	private String databaseName;
	final static Logger logger = LoggerFactory.getLogger(InterProDataLoadingServices.class);

	/**
	 * Load InterPro results list
	 * 
	 * @param databaseAccess
	 * @param map
	 * @param list
	 * @param cancel
	 * @param sequencesCounter
	 * @param auxiliaryMap 
	 */
	public InterProDataLoadingServices(Statement statement, Map<String, InterProResultsList> map,
			ConcurrentLinkedQueue<String> list, AtomicBoolean cancel, AtomicInteger sequencesCounter,
			ConcurrentHashMap<String, Integer> auxiliaryMap, String databaseName) {

		this.statement = statement;
		this.cancel = cancel;
		this.datum = sequencesCounter;
		this.interProResultsLists = map;
		this.list = list;
		this.auxiliaryMap = auxiliaryMap;
		this.databaseName = databaseName;
	}

	@Override
	public void run() {


		if(!this.cancel.get()) {
			
			try {
				
				while(!this.list.isEmpty()) {
					
					String key= this.list.poll();
					InterProResultsList resultList = this.interProResultsLists.get(key);
					
					logger.debug("key {} pool {}",key, this.list);
					
					int resultsID = InterproServices.loadInterProAnnotation(this.databaseName ,resultList.getQuery(), resultList.getQuerySequence(), 
							resultList.getMostLikelyEC(), resultList.getMostLikelyLocalization(), resultList.getName());
					
					for(InterProResult result: resultList.getResults()) {
												
						int resultID = InterproServices.loadInterProResult(this.databaseName, result.getTool(), result.geteValue(), result.getScore(), 
								result.getFamily(), result.getAccession(), result.getName(), result.getEC(), result.getGOName(),
								result.getLocalization(), result.getDatabase(), resultsID);
						
						if(result.getEntry()!=null) {
							
							int entryID = -1;
							String resultAccession = result.getEntry().getAccession();
							if(this.auxiliaryMap.containsKey(resultAccession)) {
								
								entryID = this.auxiliaryMap.get(resultAccession);
								logger.debug("Entry loaded {} {}", resultAccession, entryID);
							}
							else {								
								entryID = InterproServices.loadInterProEntry(this.databaseName,resultAccession, result.getEntry().getDescription(), 
										result.getEntry().getName(), result.getEntry().getType());
								
								this.auxiliaryMap.put(resultAccession, entryID);
							}
													
							InterproServices.loadInterProResultHasEntry(this.databaseName, resultID, entryID);
							
							for(Xref xRef : result.getEntry().getXrefs())	
								
								InterproServices.loadXrefs(this.databaseName, xRef.getCategory(), xRef.getDatabase(), xRef.getName(), xRef.getId(), entryID);
						}
						
						for(Location location: result.getLocation()) {
						
							float score = location.getScore();
							float eValue = location.getEvalue();
														
							InterproServices.loadInterProLocation(this.databaseName, location.getStart(), location.getEnd(), score, 
									location.getHmmstart(), location.getHmmend(), eValue, location.getEnvstart(),
									location.getEnvend(), location.getHmmlength(), resultID);
						}
						
						for(Model model : result.getModels()) {				
							
  							String accession = null;
  							if(this.auxiliaryMap.containsKey(model.getAccession())) {
  								
  								accession = model.getAccession();
  							}
  							else {
  								
  								accession = InterproServices.loadInterProModel(this.databaseName, model.getAccession(), model.getName(), model.getDescription());
  								
  								this.auxiliaryMap.put(accession, -1);
  								logger.debug("Model entry loaded {} ", accession);
  							}
  							
							InterproServices.loadInterProResultHasModel(this.databaseName, resultID, accession);
						}
					}
					
					InterproServices.setInterProStatus(this.databaseName, resultsID, InterproStatus.PROCESSED.toString());
					
					this.datum.incrementAndGet();
					this.setChanged();
					this.notifyObservers();
				}
			}
			catch (SQLException e) {

				e.printStackTrace();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

}
