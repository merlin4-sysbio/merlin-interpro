package pt.uminho.ceb.biosystems.merlin.interpro;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.Scanner;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.validator.routines.EmailValidator;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import es.uvigo.ei.aibench.core.operation.annotation.Cancel;
import es.uvigo.ei.aibench.core.operation.annotation.Direction;
import es.uvigo.ei.aibench.core.operation.annotation.Operation;
import es.uvigo.ei.aibench.core.operation.annotation.Port;
import es.uvigo.ei.aibench.core.operation.annotation.Progress;
import es.uvigo.ei.aibench.workbench.Workbench;
import pt.uminho.ceb.biosystems.merlin.aibench.datatypes.WorkspaceAIB;
import pt.uminho.ceb.biosystems.merlin.aibench.datatypes.annotation.AnnotationEnzymesAIB;
import pt.uminho.ceb.biosystems.merlin.aibench.gui.CustomGUI;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.MerlinUtils;
import pt.uminho.ceb.biosystems.merlin.aibench.utilities.TimeLeftProgress;
import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ebi.interpro.InterProResultsList;
import pt.uminho.ceb.biosystems.merlin.core.datatypes.WorkspaceEntity;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.InterproStatus;
import pt.uminho.ceb.biosystems.merlin.core.utilities.Enumerators.SequenceType;
import pt.uminho.ceb.biosystems.merlin.interpro.processes.InterProDomainsSearch;
import pt.uminho.ceb.biosystems.merlin.interpro.services.InterProDataLoadingServices;
import pt.uminho.ceb.biosystems.merlin.processes.annotation.AnnotationEnzymesProcesses;
import pt.uminho.ceb.biosystems.merlin.services.interpro.InterproServices;
import pt.uminho.ceb.biosystems.merlin.services.model.ModelSequenceServices;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;

/**
 * @author Oscar Dias
 *
 */
@Operation(description="operation that performs InterPro scan searches on selected sequences", name="InterPro scan")
public class HomologyInterPro implements Observer{


	private WorkspaceAIB project;
	private TimeLeftProgress progress = new TimeLeftProgress();
	private AtomicBoolean cancel;
	private Map<String,  AbstractSequence<?>> sequences;
	private double lowerThreshold;
	private double upperThreshold;
	private Set<String> genes;
	private long startTime;
	private int size;
	private String message;
	private AtomicInteger sequencesCounter;
	private int latencyWaitingPeriod;
	private boolean useManual;
	private String email;
	final static Logger logger = LoggerFactory.getLogger(HomologyInterPro.class);
	private boolean emailValid;

	@Port(direction=Direction.INPUT, name="lower threshold",description="select lower enzyme annotation threshold", defaultValue = "0", validateMethod="checkLowerThreshold", order=1)
	public void setLowerThreshold(double lowerThreshold) {

	}

	@Port(direction=Direction.INPUT, name="upper threshold",description="select upper enzyme annotation threshold", defaultValue = "0.5", validateMethod="checkUpperThreshold", order=2)
	public void setUpperThreshold(double upperThreshold) {

	}

	@Port(direction=Direction.INPUT, name="latency period",description="request latency waiting period (minimum 180 minutes)",validateMethod="checkLatencyWaitingPeriod",  advanced = true, defaultValue = "180", order=3)
	public void setLatencyWaitingPeriod(int latencyWaitingPeriod) {

		this.latencyWaitingPeriod = latencyWaitingPeriod;
	}
	
	@Port(direction=Direction.INPUT, name="check manual",description="check manually inserted entries", defaultValue = "true",  advanced = true, order=4)
	public void useManual(boolean useManual) {

		this.useManual = useManual;
	}

	/**
	 * @param project
	 * @throws Exception 
	 */
	@Port(direction=Direction.INPUT, name="workspace",description="select workspace", validateMethod="checkProject", order=6)
	public void setProject(WorkspaceAIB project) throws Exception {
		getEmail();
		if(!this.emailValid) {
			Workbench.getInstance().warn("Your email is invalid.");
			return;
		}


		if(!(this.email=="")) {
			this.project = project;

			logger.debug("Sequences {} {}", sequences.size(), sequences.keySet());

			this.sequences.keySet().retainAll(this.genes);

			logger.debug("Genes {} {}", genes.size(), genes);

			this.startTime = GregorianCalendar.getInstance().getTimeInMillis();
			this.cancel = new AtomicBoolean(false);
			this.sequencesCounter = new AtomicInteger(0);
			AtomicInteger errorCounter = new AtomicInteger(0);

			try {

				List<String> processed = InterproServices.getInterproResultsQueryByStatus(this.project.getName(), InterproStatus.PROCESSED.toString());

				logger.debug("Processed sequences {} {}", processed.size(), processed);

				this.sequences.keySet().removeAll(processed);

				logger.debug("Sequences {} {}", sequences.size(), sequences.keySet());

				if(sequences.size()>0) {

					InterProDomainsSearch interPro = new InterProDomainsSearch(this.cancel, this.sequencesCounter, errorCounter);
					interPro.addObserver(this);

					Map<String, InterProResultsList> results = null;
					this.size = this.sequences.size();
					this.message = "Searching InterPro scan";

					if(!this.cancel.get())				
						results = interPro.getInterProResults(this.sequences, (this.latencyWaitingPeriod*60*1000), this.email);

					this.startTime = GregorianCalendar.getInstance().getTimeInMillis();
					this.size = results.size();
					this.message = "Loading InterPro results";
					if(!this.cancel.get())
						this.loadInterProResults(results);

					InterproServices.deleteInterProEntries(this.project.getName(), InterproStatus.PROCESSING.toString());

					if(!this.cancel.get())	{
						if(errorCounter.get()>0) {

							Workbench.getInstance().warn("InterPro scan finished with erros. Please run again.");
						}
						else {

							MerlinUtils.updateEnzymesAnnotationView(project.getName());
							Workbench.getInstance().info("InterPro scan finished.");
						}
						
					} else {
						Workbench.getInstance().warn("InterPro scan cancelled.");
					}

				
				}
				else {

					Workbench.getInstance().info("No aminoacid sequences match the filtering criteria.");
				}
			} 
			catch (InterruptedException e) {

				Workbench.getInstance().error("InterruptedException "+e.getMessage()+" has occured.");
			}
			catch (IOException e) {

				Workbench.getInstance().error("IOException "+e.getMessage()+" has occured.");
			}
		}

	}
	
	
	/**
	 * Load the results to the database.
	 * 
	 * @param list
	 * @throws InterruptedException 
	 */
	public void loadInterProResults(Map<String, InterProResultsList> map) throws InterruptedException {

		int numberOfProcesses =  Runtime.getRuntime().availableProcessors();
		List<Thread> threads = new ArrayList<Thread>();
		List<Runnable> runnables = new ArrayList<Runnable>();
		ConcurrentLinkedQueue<String> list = new ConcurrentLinkedQueue<>(map.keySet());
		ConcurrentHashMap<String, Integer> auxiliaryMap = new ConcurrentHashMap<> ();
		
		for(int i=0; i<numberOfProcesses; i++) {

			logger.info("Starting process {} ", i);
			Runnable loadInterProData = new InterProDataLoadingServices(map, list, this.cancel, this.sequencesCounter, auxiliaryMap,this.project.getName());
			((InterProDataLoadingServices) loadInterProData).addObserver(this);
			runnables.add(loadInterProData);
			Thread thread = new Thread(loadInterProData);
			threads.add(thread);
			thread.start();
		}

		for(Thread thread :threads)			
			thread.join();
	}
	
	/**
	 * @return the progress
	 */
	@Progress
	public TimeLeftProgress getProgress() {

		return progress;
	}


	/**
	 * @param cancel the cancel to set
	 */
	@Cancel
	public void setCancel() {
		
		String[] options = new String[2];
		options[0]="yes";
		options[1]="no";

		int result=CustomGUI.stopQuestion("Cancel confirmation", "Are you sure you want to cancel the operation?", options);

		if(result==0) {
			progress.setTime(0, 0, 0);
			this.cancel.set(true);
			
			Workbench.getInstance().warn("Please hold on. Your operation is being cancelled.");
		}	
	}
	/* (non-Javadoc)
	 * @see java.util.Observer#update(java.util.Observable, java.lang.Object)
	 */
	@Override
	public void update(Observable o, Object arg) {

		this.progress.setTime(GregorianCalendar.getInstance().getTimeInMillis() - this.startTime, this.sequencesCounter.get(), this.size, this.message);
	}

	/**
	 * @param contents
	 */
	public void checkEmail(String email) {
		
			EmailValidator validator = EmailValidator.getInstance();

			if (validator.isValid(email))
				this.email = email;
			else
				throw new IllegalArgumentException("Please set a valid email address!");
	}
	
	/**
	 * @param contents
	 */
	public void checkLowerThreshold(double lowerThreshold) {

		if(lowerThreshold<0 || lowerThreshold>1)
			throw new IllegalArgumentException("the threshold should be higher than 0 and lower than 1!");

		this.lowerThreshold = lowerThreshold;
	}


	/**
	 * @param contents
	 */
	public void checkUpperThreshold(double upperThreshold) {

		if(upperThreshold<0 || upperThreshold>1)
			throw new IllegalArgumentException("the threshold should be higher than 0 and lower than 1!");

		if(upperThreshold<this.lowerThreshold)
			throw new IllegalArgumentException("the upper threshold should be higher than the lower threshold!");

		this.upperThreshold = upperThreshold;
	}

	/**
	 * @param project
	 */
	public void checkProject(WorkspaceAIB project) {

		if(project == null) {

			throw new IllegalArgumentException("no Project Selected!");
		}
		else {

			this.project = project;

			try {
				if(!ModelSequenceServices.checkGenomeSequences(project.getName(), SequenceType.PROTEIN)){//!Project.isFaaFiles(dbName, taxID)) {
					throw new IllegalArgumentException("please set the project fasta file ('.faa') to perform the InterPro similarity search.");
				}
				else {

					try {

						this.sequences = ModelSequenceServices.getGenomeFromDatabase(project.getName(), SequenceType.PROTEIN);

						if(this.sequences==null || this.sequences.isEmpty())
							throw new IllegalArgumentException("please set the project fasta files!");
					} 
					catch (Exception e) {

						throw new IllegalArgumentException("Please set the project fasta files!");
					}
				}
			}
			catch (Exception e1) {
				Workbench.getInstance().error(e1);
				e1.printStackTrace();
			}

			AnnotationEnzymesAIB enzymesAnnotation = null;

			for(WorkspaceEntity ent : project.getDatabase().getAnnotations().getEntitiesList())
				if(ent.getName().equalsIgnoreCase("enzymes"))
					enzymesAnnotation = (AnnotationEnzymesAIB) ent;

			if(enzymesAnnotation == null || enzymesAnnotation.getInitialProdItem() == null)
				throw new IllegalArgumentException("enzymes annotation view unavailable!");
			else
				this.genes = AnnotationEnzymesProcesses.getGenesInThreshold(this.lowerThreshold, this.upperThreshold, this.useManual, enzymesAnnotation);
		}
	}

	/**
	 * @param project
	 */
	public void checkLatencyWaitingPeriod(int latencyWaitingPeriod) {

		if(latencyWaitingPeriod <180)
			throw new IllegalArgumentException("The latency waiting period must be greater than 180 (zero)");

		this.latencyWaitingPeriod = latencyWaitingPeriod;
	}
	
	public void getEmail() {
		String confEmail = "";
		ArrayList <String> listLines = new ArrayList<>();
		String confPath = FileUtils.getConfFolderPath().concat("email.conf");
		File configFile = new File(confPath);
		try {
			Scanner file = new Scanner(configFile);
			while(file.hasNextLine()==true) {
				listLines.add(file.nextLine());
			}
			file.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		for (String item : listLines) {
			if(item.startsWith("email")) {
				String[] parts = item.split(":");
				confEmail = parts[1].trim();

			}

		}
		logger.debug("Email obtained from method getEmail(): " + confEmail);
		this.email = confEmail;
		
		String validation = "^[_A-Za-z0-9-\\+]+(\\.[_A-Za-z0-9-]+)*@[A-Za-z0-9-]+(\\.[A-Za-z0-9]+)*(\\.[A-Za-z]{2,})$";


		if (this.email.matches(validation)) 
			this.emailValid = true;
		else
			this.emailValid = false;
	}

}
