/**
  * Author: Dimitrios Vitsios
  * Release Date:  29 May 2013
  * Change History: -
  */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.StringTokenizer;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.biopax.paxtools.PaxtoolsMain;
import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.model.level3.Process;
import org.biopax.paxtools.io.jena.JenaIOHandler;
import org.biopax.validator.BiopaxValidatorClient;
import org.biopax.validator.BiopaxValidatorClient.RetFormat;
import org.biopax.validator.jaxb.Behavior;


public class Main {

    public static Log log = LogFactory.getLog(PaxtoolsMain.class);
    private static SimpleIOHandler io = new SimpleIOHandler();
    private static BioPAXIOHandler handler;
    private static BioPAXIOHandler handlerOut;
    private static String curDir = "";
    private static String input = "";
    private static String xmlBase = "";
    private static Model myModel;
    private static Model initialStaticModel;
    private static String pathAndOrgName = "";
    private static String pathwayName = "";
    private static String customCommentStr = "$$custom comment$$";
    private static ArrayList<String[]> customGeneIdsToNamesList;


    public static void main(String[] args) throws FileNotFoundException, IOException {

        
        curDir = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
        int lastSlashIndexInCurDir = curDir.lastIndexOf("/");
        curDir = curDir.substring(0, lastSlashIndexInCurDir);

        handlerOut = new JenaIOHandler(BioPAXLevel.L3);


        if(args.length>0){
            if(((args[0].equals("-n"))||(args[0].equals("-g"))||(args[0].equals("-i"))||(args[0].equals("-a"))||(args[0].equals("-ap"))||(args[0].equals("-an"))||(args[0].equals("-agn"))||(args[0].equals("-apgn"))||(args[0].equals("-c"))||(args[0].equals("-cp"))||(args[0].equals("-v")))&&(args.length>=2)){

                if(!(args[1].contains("/"))){
                    input = curDir+"/"+args[1];
                }
                else
                    input = args[1];

                if(!(args[0].equals("-agn")))
                    myModel = readModelL3(input);

            }
        }
        else{
            printHelp();
            System.exit(-1);
        }


        //get path_org name
        if(args[0].equals("-n") && args.length == 2)
        {
            getPathAndOrgName(input, true);
        } //getGenes
        else if(args[0].equals("-g") && args.length == 3)
        {
            getGenesFromPathway(input, args[2]);
        }
        else if(args[0].equals("-g") && args.length == 2)
        {
            getGenesFromPathway(input, "");
        }
        else if(args[0].equals("-i") && args.length == 3)
        {
            getGenesCustomIds(input, args[2], true);
        }
        else if(args[0].equals("-i") && args.length == 2)
        {
            getGenesCustomIds(input, "", true);
        } //add existence for pathway
        else if (args[0].equals("-a") && args.length == 3){
            addExistenceOfThePath(input, args[2], input);
        } //add existence for protein
        else if (args[0].equals("-ap") && args.length == 4){
            addExistenceOfTheProtein(input, args[2], args[3], true);
        } //add existence for pathway from node files
        else if (args[0].equals("-an") && args.length == 3){
            updateExistInfoForPathwayFromNodes(input, args[2]);
        } //add existence for genes/proteins from node files
        else if (args[0].equals("-agn") && args.length == 3){
            updateExistInfoForGenesFromNodes(input, args[2]);
        } //add existence for pathway and genes/proteins from node files
        else if (args[0].equals("-apgn") && args.length == 4){
            updateExistInfoForPathwayAndGenesFromNodes(input, args[2], args[3]);
        }
        //check existence of the pathway
        else if(args[0].equals("-c") && args.length == 2){
            checkExistenceOfThePath(input);
        } //check existence of a protein
        else if(args[0].equals("-cp") && args.length == 3){
            checkExistenceOfTheProtein(input, args[2]);
        } //validate the biopax pathwayNodefile
        else if(args[0].equals("-v") && ((args.length == 4)||(args.length == 3))){
            validate(args);
        }
        else if(args[0].equals("-h") && args.length == 1){
            printHelp();
        }
        else{
            printHelp();
        }



    }

    public static String getGeneSeqFromKegg(String geneId) throws IOException {

        String keggRequest = "http://rest.kegg.jp/get/"+geneId+"/aaseq";
        String geneSeq = "";

        HttpClient client = new HttpClient();
        GetMethod method = new GetMethod(keggRequest);

        // Send GET request
        int statusCode = client.executeMethod(method);

        if (statusCode != HttpStatus.SC_OK) {
          System.err.println("Method failed: " + method.getStatusLine());
        }
        InputStream rstream = null;

        // Get the response body
        rstream = method.getResponseBodyAsStream();

        // Process the response
        BufferedReader br = new BufferedReader(new InputStreamReader(rstream));
        String line;

        int lineCnt=0;

        while ((line = br.readLine()) != null) {
            if(lineCnt!=0){
                geneSeq += line+"\n";
            }
            lineCnt++;
        }
        br.close();

        return geneSeq;
   }

   public static String getGeneNameFromKegg(String geneId) throws IOException {

        String keggRequest = "http://rest.kegg.jp/get/"+geneId;
        String geneName = "";

        HttpClient client = new HttpClient();
        GetMethod method = new GetMethod(keggRequest);

        // Send GET request
        int statusCode = client.executeMethod(method);

        if (statusCode != HttpStatus.SC_OK) {
          System.err.println("Method failed: " + method.getStatusLine());
        }
        InputStream rstream = null;

        // Get the response body
        rstream = method.getResponseBodyAsStream();

        // Process the response
        BufferedReader br = new BufferedReader(new InputStreamReader(rstream));
        String line;


        while ((line = br.readLine()) != null) {

            if(line.contains("NAME")){
                StringTokenizer st = new StringTokenizer(line," ");
                st.nextToken();
                geneName += st.nextToken()+", ";
            }
            else if(line.contains("DEFINITION")){
                
                geneName += line.substring(12, line.length());
                break;
            }
        }

        br.close();

        return geneName;
   }


    public static void printHelp() throws IOException{

        System.out.println("\nInput ERROR. Use one of the following parameters format: \n" +
                            "-n [biopaxFile_input] => get pathway_organism name.\n"+
                            "-g [biopaxFile_input] ([genesList_output]) => get genes list.\n"+
                            "-i [biopaxFile_input] ([geneIds_output]) => get genes custom ids.\n"+
                            "-a [biopaxFile_input] [1/0] => add existence mark to the pathway.\n"+
                            "-ap [biopaxFile_input] [protein_id] [1/0] => add existence mark to a protein of the pathway.\n"+
                            "-an [biopaxFile_input] [pathNodeDir_path] => update pathway existence flag in a biopax file from node files information.\n"+
                            "-agn [biopaxFile_input] [geneNodesDir_path] => update proteins existence flag in a biopax file from node files information.\n"+
                            "-apgn [biopaxFile_input] [pathwayNodesDir_path] [genesNodeDir_path] => update pathway and proteins existence flag in a biopax file from node files information.\n"+
                            "-c [biopaxFile_input] => check existence mark of the pathway.\n"+
                            "-cp [biopaxFile_input] [protein_id] => check existence mark of a protein of the specified pathway.\n"+
                            "-v [biopaxFile_input] [validation_output] [xml/html] => validate the biopax file.\n"+
                            "-h -> print help.");
    }



    public static boolean getPathAndOrgName(String input, boolean printToConsole) throws FileNotFoundException{

        pathAndOrgName = "";
        boolean hasSequences = false;

        Model myModel = readModelL3(input);

        Set<Pathway> path = myModel.getObjects(Pathway.class);

        if(path.size()>1){

            int slashIndex = 0;
            if(input.contains("/"))
                slashIndex = input.lastIndexOf('/');

            pathAndOrgName = input.substring(slashIndex+1, input.length()-4);

            Iterator pathIter = path.iterator();
            while(pathIter.hasNext()){
                Pathway myPath = (Pathway)pathIter.next();
                if(myPath.getOrganism() == null)
                    continue;
                else{
                    pathwayName = myPath.getStandardName();
                }
            }

            hasSequences = false;

        } else if(path.size() == 1){


            Iterator pathIter = path.iterator();
            Pathway myPath = (Pathway)pathIter.next();
            pathwayName = myPath.getStandardName();

            
            Set<Xref> xrefs = myPath.getXref();
            Iterator xrefsIter = xrefs.iterator();

            while(xrefsIter.hasNext()){
                Xref xref = (Xref)xrefsIter.next();
                String rdfId = xref.getRDFId();
                if(rdfId.contains("Unification")){

                    pathAndOrgName += xref.getId()+"_"+xref.getDb();
                    break;
                }
            }

            hasSequences = true;
        }

        if(printToConsole)
            System.out.print(pathAndOrgName);    

        return hasSequences;
    }


    /*
     * extract genes from BioPAX pathwayNodefile and save the genes list to a FASTA formatted pathwayNodefile.
     */
    public static void getGenesFromPathway(String input, String output) throws FileNotFoundException, IOException{

        boolean hasSequences = getPathAndOrgName(input, false);

        
        ArrayList<String> protRefRdfIdsList = new ArrayList<String>();
        ArrayList<String> ecNumbersList = new ArrayList<String>();

        //get proteinReference-ECNumber mappings
        if(hasSequences){

            Set<Catalysis> catalysisElems  = myModel.getObjects(Catalysis.class);
            Iterator catalIter = catalysisElems.iterator();


            while(catalIter.hasNext()){

                Catalysis catal = (Catalysis)catalIter.next();
                System.out.println("Catalysis: "+catal.getRDFId());

                //Set<Controller> proteinsSet = catal.getController();
                for(Controller curCatalController : catal.getController()){

                    System.out.println("curCatalController: "+curCatalController.getRDFId());

                    if(curCatalController instanceof Protein){

                        Protein curProtein = (Protein)curCatalController;
                        System.out.println("curProtein (i.e. controller): "+curProtein.getRDFId());

                        for(EntityReference genericEntityReference : curProtein.getGenericEntityReferences()){

                            System.out.println("genericEntityReferences: "+genericEntityReference.getRDFId());

                            if(genericEntityReference instanceof ProteinReference){
                               String tmpProtRdfIdStr = ((ProteinReference)genericEntityReference).getRDFId();
                                protRefRdfIdsList.add(tmpProtRdfIdStr);
                                System.out.println("protRefRdfId: "+tmpProtRdfIdStr);
                            }
                        }
                    }
                    else if(curCatalController instanceof Complex){
                        Complex curComplex = (Complex)curCatalController;
                        System.out.println(">>>curComplex: "+curComplex.getRDFId());

                        Set<EntityReference> memberReferences = curComplex.getMemberReferences();
                        
                        for(EntityReference genericEntityReference: memberReferences){
                            System.out.println(">>>genericEntityReference: "+genericEntityReference.getRDFId());
                            if(genericEntityReference instanceof ProteinReference){
                               String tmpProtRdfIdStr = ((ProteinReference)genericEntityReference).getRDFId();
                                protRefRdfIdsList.add(tmpProtRdfIdStr);
                                System.out.println("protRefRdfId: "+tmpProtRdfIdStr);
                            }
                        }
                    }
                }


                Set<Process> catalProcesses = catal.getControlled();
                Iterator processIter = catalProcesses.iterator();

                while(processIter.hasNext()){
                    Process curProcess = (Process) processIter.next();
                    if(curProcess instanceof BiochemicalReaction){
                        BiochemicalReaction curReaction = (BiochemicalReaction)curProcess;
                        Set<String> curECsSet = curReaction.getECNumber();

                        Iterator ecsIter = curECsSet.iterator();
                        while(ecsIter.hasNext()){

                            String tmpCurEcStr = (String) ecsIter.next();
                            ecNumbersList.add(tmpCurEcStr);
                            System.out.println("EC: "+tmpCurEcStr+"\n");
                        }
                    }
                }
            }
        }

        
        if(output.equals("")){

            String tmpInp = input;
            tmpInp = input.substring(0, input.length()-4);
            output = tmpInp+"_FASTA.fasta";
        }
        else if(!(output.contains("/"))){
            output = curDir+"/"+output;
        }

        FileWriter foutstream = new FileWriter(output);
        BufferedWriter out = new BufferedWriter(foutstream);

        Set<ProteinReference> proteinRefs  = myModel.getObjects(ProteinReference.class);
        Iterator protIter = proteinRefs.iterator();



        Integer genesCnt = 0;

        while(protIter.hasNext())
        {
          String protAnnotStr = "";
          ProteinReference pr = (ProteinReference)protIter.next();

          String currentECNumber = "";

          for(int pId=0; pId<protRefRdfIdsList.size(); pId++){
              if(protRefRdfIdsList.get(pId).equals(pr.getRDFId())){
                currentECNumber = ecNumbersList.get(pId);
                break;
              }
          }

          Set<Xref> xrefs = pr.getXref();
          Iterator xrefIter = xrefs.iterator();


          if(hasSequences){

              genesCnt++;
              String genesCntStr = genesCnt.toString();
              String geneId = "00000";
              geneId = geneId.substring(0, (geneId.length()-genesCntStr.length()));
              geneId += genesCntStr;

              //create annotation line
              protAnnotStr += ">"+pathAndOrgName+"-"+geneId+"\t";
              while(xrefIter.hasNext()){
                  Xref xref = (Xref)xrefIter.next();
                  protAnnotStr += xref.getDb()+":";
                  protAnnotStr += xref.getId()+",";
              }
              protAnnotStr = protAnnotStr.substring(0, protAnnotStr.length()-1);
              protAnnotStr += "\t";
              protAnnotStr += pathwayName;
              protAnnotStr += "\t";
              String protNameTmp = pr.getName().toString();
              protAnnotStr += protNameTmp.substring(1, protNameTmp.length()-1);
              protAnnotStr += "\t";
              protAnnotStr += "EC:"+currentECNumber;
              protAnnotStr += "\n";
              out.write(protAnnotStr);

              
              String protSequence = pr.getSequence();
              int splitSeq;

              for(splitSeq=0; splitSeq<protSequence.length()/60; splitSeq++){
                out.write(protSequence.substring(splitSeq*60, (splitSeq+1)*60)+"\n");
              }
              out.write(protSequence.substring(splitSeq*60, protSequence.length())+"\n");
            }
            else{
                boolean includeProteinRef = false;

                while(xrefIter.hasNext()){
                    Xref xref = (Xref)xrefIter.next();
                    String rdfId = xref.getRDFId();
                    if(rdfId.contains("kegg.genes")){
                        
                          protAnnotStr = "";

                          genesCnt++;
                          String genesCntStr = genesCnt.toString();
                          String geneId = "00000";
                          geneId = geneId.substring(0, (geneId.length()-genesCntStr.length()));
                          geneId += genesCntStr;

                          //create annotation line
                          protAnnotStr += ">"+pathAndOrgName+"-"+geneId+"\t";
                          
                          protAnnotStr += xref.getDb()+":";
                          protAnnotStr += xref.getId()+",";
                         
                          protAnnotStr = protAnnotStr.substring(0, protAnnotStr.length()-1);
                          protAnnotStr += "\t";
                          protAnnotStr += pathwayName;
                          protAnnotStr += "\t";
                          String protNameTmp = pr.getName().toString();
                          protAnnotStr += protNameTmp.substring(1, protNameTmp.length()-1);
                          protAnnotStr += getGeneNameFromKegg(xref.getId())+"\n";
                          out.write(protAnnotStr);
                          out.write(getGeneSeqFromKegg(xref.getId()));
                    }
                }


            }


        }
        out.close();

        System.out.println("\nGenes list was succesfully extracted to \'"+output+"\' file.\n");
    }



    public static boolean getGenesCustomIds(String input, String output, boolean writeToFile) throws FileNotFoundException, IOException{

        customGeneIdsToNamesList = new ArrayList<String[]>();
        boolean hasSequences = getPathAndOrgName(input, false);
        FileWriter foutstream;
        BufferedWriter out = null;

        if(writeToFile){
            if(output.equals("")){

                String tmpInp = input;
                tmpInp = input.substring(0, input.length()-4);
                output = tmpInp+"_geneIds.txt";
            }
            else if(!(output.contains("/"))){
                output = curDir+"/"+output;
            }

            foutstream = new FileWriter(output);
            out = new BufferedWriter(foutstream);
        }

        Set<ProteinReference> proteinRefs  = myModel.getObjects(ProteinReference.class);
        Iterator protIter = proteinRefs.iterator();
        int counterRef = 0;


        Integer genesCnt = 0;

        while(protIter.hasNext())
        {

          ProteinReference pr = (ProteinReference)protIter.next();

          Set<Xref> xrefs = pr.getXref();
          Iterator xrefIter = xrefs.iterator();


          String protName = "";

          if(hasSequences){

              genesCnt++;
              String genesCntStr = genesCnt.toString();
              String geneId = "00000";
              geneId = geneId.substring(0, (geneId.length()-genesCntStr.length()));
              geneId += genesCntStr;
              geneId = pathAndOrgName+"-"+geneId;

              while(xrefIter.hasNext()){
                  Xref xref = (Xref)xrefIter.next();
                  protName += xref.getDb()+":";
                  protName += xref.getId()+",";
              }
              protName = protName.substring(0, protName.length()-1);

              
              String[] tmpStringArray = new String[2];
              tmpStringArray[0] = geneId;
              tmpStringArray[1] = protName;

              customGeneIdsToNamesList.add(tmpStringArray);

            }
            else{

                while(xrefIter.hasNext()){
                    Xref xref = (Xref)xrefIter.next();
                    String rdfId = xref.getRDFId();
                    if(rdfId.contains("kegg.genes")){
                        
                          protName = "";

                          genesCnt++;
                          String genesCntStr = genesCnt.toString();
                          String geneId = "00000";
                          geneId = geneId.substring(0, (geneId.length()-genesCntStr.length()));
                          geneId += genesCntStr;
                          geneId = pathAndOrgName+"-"+geneId;

                          protName += xref.getDb()+":";
                          protName += xref.getId();


                          String[] tmpStringArray = new String[2];
                          tmpStringArray[0] = geneId;
                          tmpStringArray[1] = protName;

                          customGeneIdsToNamesList.add(tmpStringArray);
                          
                    }
                }

            }


        }

       

        if(writeToFile){

            for(int i=0; i<customGeneIdsToNamesList.size(); i++){
                String[] tmpArr = customGeneIdsToNamesList.get(i);
                out.write(tmpArr[0]+"\n");
            }
        
            
            out.close();

            System.out.println("\nGenes IDs were succesfully extracted to \'"+output+"\' file.\n");
        }

        return hasSequences;
    }



    public static void addExistenceOfThePath(String inputPath, String existSign, String outputPath) throws IOException{

        boolean inputEqualsOutput = inputPath.equals(outputPath);

        Set<Pathway> pwSet = myModel.getObjects(Pathway.class);
        Iterator itPath = pwSet.iterator();
        Pathway pw = (Pathway)itPath.next();


        if(Integer.parseInt(existSign) == 1){


            boolean addCommentFlag = true;
            Set<String> coms = pw.getComment();
            Iterator itComs = coms.iterator();

            while(itComs.hasNext()){

                String com = (String)itComs.next();
                if(com.equals(customCommentStr+":present")){
                    if(inputEqualsOutput)
                        System.out.println("\nThe pathway is already set as 'present'!\n");
                    addCommentFlag = false;
                    break;
                }
                else if(com.equals(customCommentStr+":absent")){
                    pw.removeComment(customCommentStr+":absent");
                    break;
                }

            }

            if(addCommentFlag){
                pw.addComment(customCommentStr+":present");
                if(inputEqualsOutput)
                    System.out.println("\nThe pathway was set as 'present'.\n");
            }

        }
        else if(Integer.parseInt(existSign) == 0){

            boolean addCommentFlag = true;
            Set<String> coms = pw.getComment();
            Iterator itComs = coms.iterator();

            while(itComs.hasNext()){

                String com = (String)itComs.next();
                if(com.equals(customCommentStr+":absent")){
                    if(inputEqualsOutput)
                        System.out.println("\nThe pathway is already set as 'absent'!\n");
                    addCommentFlag = false;
                    break;
                }
                else if(com.equals(customCommentStr+":present")){
                    pw.removeComment(customCommentStr+":present");
                    break;
                }

            }

            if(addCommentFlag){
                pw.addComment(customCommentStr+":absent");
                if(inputEqualsOutput)
                    System.out.println("\nThe pathway was set as 'absent'.\n");
            }
        }
        else{
            System.out.println("Input ERROR. Correct format:\n"+
                    "-a [biopaxFile_input] [1/0]\n");
        }

        OutputStream outputToOwl = new FileOutputStream(outputPath);
        handlerOut.convertToOWL(myModel, outputToOwl);
    }


    public static void checkExistenceOfThePath(String inputPath){

        Set<Pathway> pwSet = myModel.getObjects(Pathway.class);
        Iterator itPath = pwSet.iterator();
        Pathway pw = (Pathway)itPath.next();
        

        Set<String> pwComs = pw.getComment();
        Iterator itPwComs = pwComs.iterator();

        boolean existAttrSet = false;

        int cnt = 0;
        while(itPwComs.hasNext()){

            String tmpCom = (String)itPwComs.next();

            if(tmpCom.equals(customCommentStr+":present")){
                System.out.println("\n[1]: The pathway is present!\n");
                existAttrSet = true;
                break;
            }
            else if(tmpCom.equals(customCommentStr+":absent")){
                System.out.println("\n[0]: The pathway is absent!\n");
                existAttrSet = false;
                break;
            }
        }

        if(!existAttrSet)
            System.out.println("\nNo existence property (presence/absence) is set for the specified pathway.");
    }

    public static void checkExistenceOfTheProtein(String inputPath, String protCustomId) throws FileNotFoundException, IOException{

        boolean hasSequences = getGenesCustomIds(inputPath,"", false);

        Set<ProteinReference> prSet = myModel.getObjects(ProteinReference.class);
        Iterator itProt = prSet.iterator();

        
        boolean validProteinName = false;
        boolean foundAtLeastOne = false;

        while(itProt.hasNext()){

            String protNameTmp = "";
            boolean foundProtein = false;
            ProteinReference pr = (ProteinReference)itProt.next();

            if(hasSequences){
                Set<Xref> xrefs = pr.getXref();
                Iterator xrefIter = xrefs.iterator();


                while(xrefIter.hasNext()){
                      Xref xref = (Xref)xrefIter.next();
                      protNameTmp += xref.getDb()+":";
                      protNameTmp += xref.getId()+",";
                  }
                  protNameTmp = protNameTmp.substring(0, protNameTmp.length()-1);


                for(int listrow=0; listrow<customGeneIdsToNamesList.size(); listrow++){
                    String[] tmpRow = customGeneIdsToNamesList.get(listrow);
                    if(tmpRow[0].equals(protCustomId) && tmpRow[1].equals(protNameTmp)){
                        foundProtein = true;
                        break;
                    }
                }
            }
            else{

                Set<Xref> xrefs = pr.getXref();
                Iterator xrefIter = xrefs.iterator();

                while(xrefIter.hasNext()){
                    Xref xref = (Xref)xrefIter.next();
                    String rdfId = xref.getRDFId();
                    if(rdfId.contains("kegg.genes")){

                          protNameTmp += xref.getDb()+":";
                          protNameTmp += xref.getId();


                          for(int listrow=0; listrow<customGeneIdsToNamesList.size(); listrow++){
                            String[] tmpRow = customGeneIdsToNamesList.get(listrow);
                            if(tmpRow[0].equals(protCustomId) && tmpRow[1].equals(protNameTmp)){
                                foundProtein = true;
                                break;
                            }
                        }

                    }
                }

            }

            if(foundProtein){

                validProteinName = true;
                Set<String> coms = pr.getComment();
                Iterator itCom = coms.iterator();

                while(itCom.hasNext()){

                    String tmpCom = (String)itCom.next();

                    if(tmpCom.equals(customCommentStr+":present")){
                        System.out.println("\n[1]: The protein is present!\n");
                        foundAtLeastOne = true;
                        break;
                    }
                    else if(tmpCom.equals(customCommentStr+":absent")){
                        System.out.println("\n[0]: The protein is absent!\n");
                        foundAtLeastOne = true;
                        break;
                    }
                }
            }
        }

        if(!validProteinName)
            System.out.println("\nNo such protein at the specified pathway.");
        else if(!foundAtLeastOne)
            System.out.println("\nNo existence attribute (presence/absence) is set for the specified protein.");

    }


    public static void addExistenceOfTheProtein(String inputPath, String protCustomId, String existSign, boolean singleOperatingCommand) throws IOException{

        
        boolean hasSequences;


        if(singleOperatingCommand)
            hasSequences = getGenesCustomIds(inputPath,"", false);
        else
            hasSequences = getPathAndOrgName(inputPath,false);

        Set<ProteinReference> prSet = myModel.getObjects(ProteinReference.class);
        Iterator itProt = prSet.iterator();

        
        boolean foundAtLeastOne = false;

        while(itProt.hasNext()){

            boolean foundProtein = false;
            String protNameTmp = "";
            ProteinReference pr = (ProteinReference)itProt.next();

            if(hasSequences){

                Set<Xref> xrefs = pr.getXref();
                Iterator xrefIter = xrefs.iterator();

                while(xrefIter.hasNext()){
                      Xref xref = (Xref)xrefIter.next();
                      protNameTmp += xref.getDb()+":";
                      protNameTmp += xref.getId()+",";
                }
                protNameTmp = protNameTmp.substring(0, protNameTmp.length()-1);


                for(int listrow=0; listrow<customGeneIdsToNamesList.size(); listrow++){
                    String[] tmpRow = customGeneIdsToNamesList.get(listrow);
                    if(tmpRow[0].equals(protCustomId) && tmpRow[1].equals(protNameTmp)){
                        foundProtein = true;
                        break;
                    }
                }
            }
            else{

                Set<Xref> xrefs = pr.getXref();
                Iterator xrefIter = xrefs.iterator();

                while(xrefIter.hasNext()){
                    Xref xref = (Xref)xrefIter.next();
                    String rdfId = xref.getRDFId();
                    if(rdfId.contains("kegg.genes")){

                          protNameTmp += xref.getDb()+":";
                          protNameTmp += xref.getId();


                          for(int listrow=0; listrow<customGeneIdsToNamesList.size(); listrow++){
                            String[] tmpRow = customGeneIdsToNamesList.get(listrow);
                            if(tmpRow[0].equals(protCustomId) && tmpRow[1].equals(protNameTmp)){
                                foundProtein = true;
                                break;
                            }
                        }

                    }
                }

            }

           

            if(foundProtein){

                if(Integer.parseInt(existSign) == 1){


                    boolean addCommentFlag = true;
                    Set<String> coms = pr.getComment();
                    Iterator itComs = coms.iterator();

                    while(itComs.hasNext()){

                        String com = (String)itComs.next();


                        if(com.equals(customCommentStr+":present")){
                            System.out.println("\nThe protein is already set as 'present'!");
                            addCommentFlag = false;
                            break;
                        }
                        else if(com.equals(customCommentStr+":absent")){
                            pr.removeComment(customCommentStr+":absent");
                            break;
                        }

                    }

                    if(addCommentFlag){
                        pr.addComment(customCommentStr+":present");
                        System.out.println("\nThe protein was set as 'present'.");
                    }
                    

                }
                else if(Integer.parseInt(existSign) == 0){

                    boolean addCommentFlag = true;
                    Set<String> coms = pr.getComment();
                    Iterator itComs = coms.iterator();

                    while(itComs.hasNext()){

                        String com = (String)itComs.next();
                        
                        if(com.equals(customCommentStr+":absent")){
                            System.out.println("\nThe protein is already set as 'absent'!");
                            addCommentFlag = false;
                            break;
                        }
                        else if(com.equals(customCommentStr+":present")){
                            pr.removeComment(customCommentStr+":present");
                            break;
                        }

                    }

                    if(addCommentFlag){
                        pr.addComment(customCommentStr+":absent");
                        System.out.println("\nThe protein was set as 'absent'.");
                    } 
                }
                else{
                    System.out.println("Input ERROR. Correct format:\n"+
                            "-ap [biopaxFile_input] [protein_id] [1/0]");
                }


                foundAtLeastOne = true;
                break;
            }
        }

        if(!foundAtLeastOne)
            System.out.println("\nNo such protein in the specified pathway!");

        if(singleOperatingCommand){
            OutputStream outputToOwl = new FileOutputStream(inputPath);
            handlerOut.convertToOWL(myModel, outputToOwl);
        }
    }


    public static void updateExistInfoForPathwayFromNodes(String biopaxFilePath, String pathNodesDirPath) throws FileNotFoundException{


        pathAndOrgName = "";
        getPathAndOrgName(biopaxFilePath, false);
        System.out.println("-Pathway name: "+pathAndOrgName+"\n");

        String outputDir = curDir + "/output";
        String newBiopaxListDir = outputDir+"/"+pathAndOrgName;
        boolean createOuputDir = new File(outputDir).mkdir();
        boolean createNewBiopaxListDir = new File(newBiopaxListDir).mkdir();



        File nodeF = new File(pathNodesDirPath);
        File[] nodeFiles = nodeF.listFiles();
        for (File pathwayNodefile : nodeFiles){

            try{

                System.out.println("> Checking nodefile: "+pathwayNodefile.getName());
                FileInputStream pwyFstream = new FileInputStream(pathNodesDirPath+"/"+pathwayNodefile.getName());

                DataInputStream pwyIn = new DataInputStream(pwyFstream);
                BufferedReader pwyBr = new BufferedReader(new InputStreamReader(pwyIn));

                String strLine = "";
                String pathNameFromNode = "";

                strLine = pwyBr.readLine();

                boolean addedAtLeastByOneNode = false;

                while (!((strLine = pwyBr.readLine()).contains("Total number of analysed families"))){
                    StringTokenizer st = new StringTokenizer(strLine);
                    pathNameFromNode = st.nextToken();


                    if(pathAndOrgName.equals(pathNameFromNode)){
                        addedAtLeastByOneNode = true;
                        break;
                    }

                }

                String outputPath = newBiopaxListDir;
                outputPath += "/"+(pathwayNodefile.getName()).substring(0, (pathwayNodefile.getName()).length()-5)+"_";
                outputPath += pathAndOrgName+".owl";
                if(addedAtLeastByOneNode){
                    addExistenceOfThePath(biopaxFilePath, "1", outputPath);
                    System.out.println("pathway 'presence' property was set.\n");
                }
                else {
                    addExistenceOfThePath(biopaxFilePath, "0", outputPath);
                    System.out.println("pathway 'absence' property was set.\n");
                }


                pwyIn.close();
            } catch (Exception e){
                System.err.println("Error: " + e.getMessage());
            }

        }

    }


    public static void updateExistInfoForPathwayAndGenesFromNodes(String biopaxFilePath, String pathNodesDirPath, String geneNodesDirPath) throws FileNotFoundException, IOException{

        pathAndOrgName = "";
        getPathAndOrgName(biopaxFilePath, false);
        System.out.println("-Pathway name: "+pathAndOrgName+"\n");

        String outputDir = curDir + "/output";
        String newBiopaxListDir = outputDir+"/"+pathAndOrgName;
        boolean createOuputDir = new File(outputDir).mkdir();
        boolean createNewBiopaxListDir = new File(newBiopaxListDir).mkdir();


        File nodeF = new File(geneNodesDirPath);
        File[] nodeFiles = nodeF.listFiles();
        for (File file : nodeFiles){

            ArrayList<String> presentGenes = new ArrayList<String>();

            String outputPath = newBiopaxListDir;
            outputPath += "/"+(file.getName()).substring(0, (file.getName()).length()-5)+"_";
            getPathAndOrgName(biopaxFilePath, false);
            outputPath += pathAndOrgName+".owl";

            
            try{
                File f1 = new File(biopaxFilePath);
                File f2 = new File(outputPath);
                InputStream in = new FileInputStream(f1);

                OutputStream out = new FileOutputStream(f2);

                byte[] buf = new byte[1024];
                int len;
                while ((len = in.read(buf)) > 0){
                    out.write(buf, 0, len);
                }

                in.close();
                out.close();
                System.out.println("File copied.");
            }
            catch(FileNotFoundException ex){

                System.out.println(ex.getMessage() + " in the specified directory.");
                System.exit(0);
            }
            catch(IOException e){
                System.out.println("-apgn breakpoint 1:"+ e.getMessage());
            }

            myModel = readModelL3(biopaxFilePath);

            getGenesCustomIds(biopaxFilePath,"", false);

            try{

                System.out.println("> Checking nodefile: "+file.getName());
                FileInputStream fstream = new FileInputStream(geneNodesDirPath+"/"+file.getName());

                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));

                String strLine = "";
                String geneNameFromNode = "";

                strLine = br.readLine();

                boolean addedAtLeastByOneNode = false;

                while (!((strLine = br.readLine()).contains("Total number of analysed families"))){
                    StringTokenizer st = new StringTokenizer(strLine);
                    geneNameFromNode = st.nextToken();

                    addExistenceOfTheProtein(biopaxFilePath, geneNameFromNode, "1", false);
                    presentGenes.add(geneNameFromNode);
                }

                in.close();
                fstream.close();


                for(int listrow=0; listrow<customGeneIdsToNamesList.size(); listrow++){
                    String[] tmpRow = customGeneIdsToNamesList.get(listrow);
                    String tmpCustomId = tmpRow[0];

                    boolean isPresent = false;

                    for(int presGenId = 0; presGenId<presentGenes.size(); presGenId++){
                        if((presentGenes.get(presGenId)).equals(tmpCustomId)){
                            isPresent = true;
                            break;
                        }
                    }

                    if(!isPresent){
                        addExistenceOfTheProtein(outputPath, tmpRow[0], "0", false);
                    }

                }
                
                

            } catch (Exception e){
                System.err.println("-apgn breakpoint 2, Error: " + e.getMessage());
            }

            try{
                FileInputStream pwyFstream = new FileInputStream(pathNodesDirPath+"/"+file.getName());

                DataInputStream pwyIn = new DataInputStream(pwyFstream);
                BufferedReader pwyBr = new BufferedReader(new InputStreamReader(pwyIn));

                String strLine = "";
                String pathNameFromNode = "";

                strLine = pwyBr.readLine();

                boolean addedAtLeastByOneNode = false;

                while (!((strLine = pwyBr.readLine()).contains("Total number of analysed families"))){
                    StringTokenizer st = new StringTokenizer(strLine);
                    pathNameFromNode = st.nextToken();


                    if(pathAndOrgName.equals(pathNameFromNode)){
                        addedAtLeastByOneNode = true;
                        break;
                    }

                }

                if(addedAtLeastByOneNode){
                    addExistenceOfThePath(biopaxFilePath, "1", outputPath);
                    System.out.println("pathway 'presence' property was set.\n");
                }
                else  {
                    addExistenceOfThePath(biopaxFilePath, "0", outputPath);
                    System.out.println("pathway 'absence' property was set.\n");
                }


                pwyIn.close();
            }  catch (Exception e){
                System.err.println("-apgn breakpoint 3: Error: " + e.getMessage());
            }

        }

    }

    public static void updateExistInfoForGenesFromNodes(String biopaxFilePath, String geneNodesDirPath) throws FileNotFoundException, IOException {

        pathAndOrgName = "";
        getPathAndOrgName(biopaxFilePath, false);
        System.out.println("-Pathway name: "+pathAndOrgName+"\n");

        String outputDir = curDir + "/output";
        String newBiopaxListDir = outputDir+"/"+pathAndOrgName;
        boolean createOuputDir = new File(outputDir).mkdir();
        boolean createNewBiopaxListDir = new File(newBiopaxListDir).mkdir();


        File nodeF = new File(geneNodesDirPath);
        File[] nodeFiles = nodeF.listFiles();
        for (File file : nodeFiles){

            ArrayList<String> presentGenes = new ArrayList<String>();

            String outputPath = newBiopaxListDir;
            outputPath += "/"+(file.getName()).substring(0, (file.getName()).length()-5)+"_";
            outputPath += pathAndOrgName+".owl";

            
            try{
                File f1 = new File(biopaxFilePath);
                File f2 = new File(outputPath);
                InputStream in = new FileInputStream(f1);

                OutputStream out = new FileOutputStream(f2);

                byte[] buf = new byte[1024];
                int len;
                while ((len = in.read(buf)) > 0){
                    out.write(buf, 0, len);
                }
                
                in.close();
                out.close();
                System.out.println("File copied.");
            }
            catch(FileNotFoundException ex){
                
                System.out.println(ex.getMessage() + " in the specified directory.");
                System.exit(0);
            }
            catch(IOException e){
                System.out.println(e.getMessage());  
            }

            myModel = readModelL3(biopaxFilePath);

           
            getGenesCustomIds(biopaxFilePath,"", false);

            
            try{

                System.out.println("> Checking nodefile: "+file.getName());
                FileInputStream fstream = new FileInputStream(geneNodesDirPath+"/"+file.getName());

                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));

                String strLine = "";
                String geneNameFromNode = "";

                strLine = br.readLine();

                boolean addedAtLeastByOneNode = false;

                while (!((strLine = br.readLine()).contains("Total number of analysed families"))){
                    StringTokenizer st = new StringTokenizer(strLine);
                    geneNameFromNode = st.nextToken();


                    addExistenceOfTheProtein(biopaxFilePath, geneNameFromNode, "1", false);
                    
                    presentGenes.add(geneNameFromNode);
                }

                in.close();
                fstream.close();

                for(int listrow=0; listrow<customGeneIdsToNamesList.size(); listrow++){
                    String[] tmpRow = customGeneIdsToNamesList.get(listrow);
                    String tmpCustomId = tmpRow[0];

                    boolean isPresent = false;

                    for(int presGenId = 0; presGenId<presentGenes.size(); presGenId++){
                        if((presentGenes.get(presGenId)).equals(tmpCustomId)){
                            isPresent = true;
                            break;
                        }
                    }

                    if(!isPresent){
                        
                        addExistenceOfTheProtein(biopaxFilePath, tmpRow[0], "0", false);
                    }

                }

               

            } catch (Exception e){
                System.err.println("Error: " + e.getMessage());
            }

            OutputStream outputToOwl = new FileOutputStream(outputPath);
            handlerOut.convertToOWL(myModel, outputToOwl);
            outputToOwl.close();
        }




    }


    public static Model readModelL3(String inputFileStr) throws FileNotFoundException{

        handler = new JenaIOHandler(BioPAXLevel.L3);

        InputStream inputStreamFromFile = new FileInputStream(inputFileStr);

        Model model = handler.convertFromOWL(inputStreamFromFile);

        return model;
    }




    /**
     * Checks files by the online BioPAX validator
     * using the validator client.
     *
     * @see <a href="http://www.biopax.org/biopax-validator/ws.html">BioPAX Validator Webservice</a>
     *
     * @param argv
     * @throws IOException
     */
    public static void validate(String[] argv) throws IOException
    {
        String input = argv[1];
        String output = argv[2];
        // default options
        RetFormat outf = RetFormat.XML;
        boolean fix = false;
        boolean normalize = false;
        Integer maxErrs = null;
        Behavior level = null; //will report both errors and warnings
        // match optional args
		for (int i = 3; i < argv.length; i++) {
			if ("html".equalsIgnoreCase(argv[i])) {
				outf = RetFormat.HTML;
			} else if ("xml".equalsIgnoreCase(argv[i])) {
				outf = RetFormat.XML;
			} else if ("biopax".equalsIgnoreCase(argv[i])) {
				outf = RetFormat.OWL;
			} else if ("normalize".equalsIgnoreCase(argv[i])) {
				normalize = true;
			} else if ("auto-fix".equalsIgnoreCase(argv[i])) {
				fix = true;
			} else if ("only-errors".equalsIgnoreCase(argv[i])) {
				level = Behavior.ERROR;
			} else if ((argv[i]).toLowerCase().startsWith("maxerrors=")) {
				String num = argv[i].substring(10);
				maxErrs = Integer.valueOf(num);
			}
		}

        Collection<File> files = new HashSet<File>();
        File fileOrDir = new File(input);
        if (!fileOrDir.canRead()) {
            System.out.println("Cannot read " + input);
            System.exit(-1);
        }

        // collect files
        if (fileOrDir.isDirectory()) {
            // validate all the OWL files pwyIn the folder
            FilenameFilter filter = new FilenameFilter() {
                public boolean accept(File dir, String name) {
                    return (name.endsWith(".owl"));
                }
            };
            for (String s : fileOrDir.list(filter)) {
                files.add(new File(fileOrDir.getCanonicalPath()
                        + File.separator + s));
            }
        } else {
            files.add(fileOrDir);
        }

        // upload and validate using the default URL:
        // http://www.biopax.org/biopax-validator/check.html
        OutputStream os = new FileOutputStream(output);
        try {
            if (!files.isEmpty()) {
                BiopaxValidatorClient val =
                        new BiopaxValidatorClient("http://www.biopax.org/biopax-validator/check.html");

                ByteArrayOutputStream baos = new ByteArrayOutputStream();
                val.validate(fix, normalize, outf, level, maxErrs, null, files.toArray(new File[]{}), os);
                //System.out.println(os.toString());
                // do not auto-fix, nor normalize, nor filter errors, etc..
//                val.validate(fix, normalize, outf, level, maxErrs, null,
// files.toArray(new File[]{}), os); //Todo ED>>Igor>>Does not compile on my
// system. Prob. dep issue?
            }
        } catch (Exception ex) {
            // fall-back: not using the remote validator; trying to read files
            String msg = "Faild to Validate Using the Remote Service.\n " +
                    "Now Trying To Read Each File and Build The Model\n" +
                    "Watch Log Messages...\n";
            System.err.println(msg);
            os.write(msg.getBytes());

            for (File f : files) {
                Model m = null;
                msg = "";
                try {
                    m = io.convertFromOWL(new FileInputStream(f));
                    msg = "Model that contains "
                            + m.getObjects().size()
                            + " elements is created (check the log messages)\n";
                    os.write(msg.getBytes());
                } catch (Exception e) {
                    msg = "Error during validation" + e + "\n";
                    os.write(msg.getBytes());
                    e.printStackTrace();
                    log.error(msg);
                }
                os.flush();
            }
        }
    }





}

