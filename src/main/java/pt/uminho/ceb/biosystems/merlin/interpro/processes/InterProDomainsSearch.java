package pt.uminho.ceb.biosystems.merlin.interpro.processes;

import java.io.IOException;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.EbiAPI;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.interpro.InterProResultsList;

/**
 * @author Oscar Dias
 *
 */
public class InterProDomainsSearch extends Observable implements Observer {

	private AtomicBoolean cancel;
	private AtomicInteger errorCounter;
	private AtomicInteger sequencesCounter;
	final static Logger logger = LoggerFactory.getLogger(InterProDomainsSearch.class);
	
	/**
	 * Constructor for processing InterPro scan.
	 * 
	 * @param databaseAccess
	 * @param cancel
	 * @param sequencesCounter
	 * @param errorCounter 
	 */
	public InterProDomainsSearch(AtomicBoolean cancel, AtomicInteger sequencesCounter, AtomicInteger errorCounter) {
		
		this.cancel = cancel;
		this.sequencesCounter = sequencesCounter;
		this.errorCounter = errorCounter;
	}
	
	/**
	 * Get InterPro results map.
	 * @param waitingPeriod 
	 * 
	 * @return
	 * @throws InterruptedException 
	 * @throws IOException 
	 * @throws Exception
	 */
	public Map<String, InterProResultsList> getInterProResults(Map<String, AbstractSequence<?>> genome, long waitingPeriod, String email) throws InterruptedException, IOException  {
		
		EbiAPI ebiAPI = new EbiAPI();
		ebiAPI.addObserver(this);
		Map<String, InterProResultsList> interPro = ebiAPI.getInterProAnnotations(genome, this.errorCounter, this.cancel, waitingPeriod, this.sequencesCounter, email);
		
		return interPro;
	}
	
	@Override
	public void update(java.util.Observable o, Object arg) {
		
		logger.debug("sequence updated {}", this.sequencesCounter.get());
		
		setChanged();
		notifyObservers();
	}

}
