package crawler.com;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;


import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

import com.gargoylesoftware.htmlunit.FailingHttpStatusCodeException;
import com.gargoylesoftware.htmlunit.WebClient;
import com.gargoylesoftware.htmlunit.html.HtmlAnchor;
import com.gargoylesoftware.htmlunit.html.HtmlPage;

public class DBASSCrawler {
	private static final String URL = "http://www.dbass.org.uk/DBASS5/viewlist.aspx";
	private static final String BASEURL = "http://www.dbass.org.uk/DBASS5";
	private static ArrayList<String> links= new ArrayList<String>();
	private static ArrayList<String> crypticSpliceSites= new ArrayList<String>();
	private static ArrayList<String> neighboringSpliceSites= new ArrayList<String>();
	private static HashMap<String, String> authenticAndCrypticSpliSite= new HashMap<String, String>();
	public static void getPageLinks(WebClient webClient,HtmlPage page) {
		try  {
			String pageAsXml = page.asXml();
			//HtmlAnchor anchor = page.getAnchorByName("anchor_name");
			Document document = Jsoup.parse(pageAsXml);
			Elements linksOnPage = document.select("a[title=View details for this record]");
				                
            for (Element link : linksOnPage) {
            	if (!links.contains(link.attr("href").substring(1))) {
            		links.add(link.attr("href").substring(1));
                    //System.out.println(link.attr("href"));
                }
            }
            HtmlAnchor nextPage = page.getAnchorByText("[Next Page]");               
            page = nextPage.click();
            webClient.waitForBackgroundJavaScript(1000);
            document = Jsoup.parse(page.asXml());
            Elements currentPage = document.select("span[id=PageBody_lblCurrentPage]");
            Elements lastPage = document.select("span[id=PageBody_lblTotalPages]");
            if(Integer.parseInt(currentPage.first().text()) < Integer.parseInt(lastPage.first().text())){
            	getPageLinks(webClient,page);
            }
            webClient.close();
	    } catch (FailingHttpStatusCodeException e) {
			
			e.printStackTrace();
		} catch (MalformedURLException e) {
			
			e.printStackTrace();
		} catch (IOException e) {
			
			e.printStackTrace();
		}
	}
	public static  void get9mers(){
		Document document;
		String nucleotideSequence;
		int distBetweenAuthAbbr=Integer.MIN_VALUE;
		//String link = "/viewsplicesite.aspx?id=388";
		for (String link : links) {
			//System.out.println(link);
			nucleotideSequence="";
			//distBetweenAuthAbbr =0;
			try{
				String childURL = BASEURL+link;
				document= Jsoup.connect(childURL).get();
				Element table = document.select("table.Dbass").first();
				Elements trs = table.getElementsByTag("tr");
				for(Element tr : trs){
					boolean flag =false;
					for (Element td : tr.children()){
						if(td.text().equalsIgnoreCase("Distance between Authentic and Aberrant 5' Splice Site (nt)")){
							//System.out.println(tr.getElementsByClass("center").text());
							if(!tr.getElementsByClass("center").text().equals("-")){
								distBetweenAuthAbbr = Integer.parseInt(tr.getElementsByClass("center").text());
							}
							flag = true;
							break;
						}
					}
					if(flag){
						break;
					}
				}
				Elements divPSequence = document.select("div[id=PageBody_pnlSequence]");
				for (Element element : divPSequence) {
					for(Element span:element.getElementsByTag("span")){
						if(span.hasClass("exon")){
							nucleotideSequence=nucleotideSequence.concat(span.getElementsByClass("exon").text());
							//System.out.println(span.getElementsByClass("exon").text());
						}
						else if(span.hasClass("intron")){
							nucleotideSequence=nucleotideSequence.concat(span.getElementsByClass("intron").text());
							//System.out.println(span.getElementsByClass("intron").text());
							
						}
						else if(span.hasClass("marker")){
							nucleotideSequence=nucleotideSequence.concat(span.getElementsByClass("marker").text());
						}
						else if(span.hasClass("pseudoexon")){
							nucleotideSequence=nucleotideSequence.concat(span.getElementsByClass("pseudoexon").text());
						}
					}		
				}		
				getCrypticAuthentic9Mer(nucleotideSequence.toUpperCase(),distBetweenAuthAbbr,link);
				
			}catch(IOException ex)
			{
				ex.printStackTrace();
			}		
		}
	}
	public static void getCrypticAuthentic9Mer(String nucleotideSequence,int distBetweenAuthAbbr,String link){	
		String authentic9Mer="ZZZZZZZZZ",cryptic9Mer = "ZZZZZZZZZ";
		String nucleotideAfterMutation="",nucleotideBeforeMutation="";
		HashMap<String, String> nucleotideAuthCrpt = checkDeletionMutation(nucleotideSequence);
		for ( Entry<String, String> entry : nucleotideAuthCrpt.entrySet()) {
			nucleotideBeforeMutation= entry.getKey();
			nucleotideAfterMutation = entry.getValue();
		}
		int crypticSpliceSite = nucleotideAfterMutation.indexOf("/GT");
		int authenticSpliceSite = Integer.MAX_VALUE;
		
		if(crypticSpliceSite>0){
			/*System.out.println(nucleotideSequence);
			System.out.println(nucleotideBeforeMutation);
			System.out.println(nucleotideAfterMutation); */
			
			if(distBetweenAuthAbbr!=Integer.MIN_VALUE){
				if(distBetweenAuthAbbr < 0){
					authenticSpliceSite = crypticSpliceSite+1-distBetweenAuthAbbr;
				}
				else{
					authenticSpliceSite = crypticSpliceSite-distBetweenAuthAbbr;
				}
				if(authenticSpliceSite>=3 && authenticSpliceSite<(nucleotideBeforeMutation.length()-6)){	
					boolean check = nucleotideBeforeMutation.substring(authenticSpliceSite-3, authenticSpliceSite+6).contains("/");
					if(check){
						String temp = nucleotideBeforeMutation.replace("/", "");
						authentic9Mer = temp.substring(authenticSpliceSite-4, authenticSpliceSite+5);
					}
					else
					{
						authentic9Mer = nucleotideBeforeMutation.substring(authenticSpliceSite-3, authenticSpliceSite+6);
					}
					//System.out.println("Authrntic 9-mer: "+authentic9Mer);
				}
				//System.out.println(crypticSpliceSite+"*****"+authenticSpliceSite);
			}
			cryptic9Mer = nucleotideAfterMutation.substring(crypticSpliceSite-3, crypticSpliceSite).concat(nucleotideAfterMutation.substring(crypticSpliceSite+1, crypticSpliceSite+7));
			//System.out.println("Cyptic 9-mer: "+cryptic9Mer);
			crypticSpliceSites.add(cryptic9Mer);
			authenticAndCrypticSpliSite.put(authentic9Mer, cryptic9Mer);
			getNeighborData(nucleotideAfterMutation,link,cryptic9Mer,authentic9Mer);
		}
		else
		{
			//System.out.println("Sites having issues with Cryptic splicesite: "+link);
		}
	}
	public static HashMap<String,String> checkDeletionMutation(String nucleotideSequence){
		HashMap<String,String> nucleotideAuthCrpt = new HashMap<String,String>();
		String stringTobeReplaced ="";
		String nucleotideAfterMutation =nucleotideSequence,nucleotideBeforeMutation=nucleotideSequence;
		if(nucleotideSequence.contains(">")){
			//System.out.println("There has been a mutation of type: Substitution");
			stringTobeReplaced = nucleotideSequence.substring(nucleotideSequence.indexOf('(')-1, nucleotideSequence.indexOf(')')+1);
			if(stringTobeReplaced.equals("/()")){
				nucleotideSequence=nucleotideSequence.replace(stringTobeReplaced, "");
			}
			nucleotideAfterMutation =nucleotideSequence;nucleotideBeforeMutation=nucleotideSequence;
			stringTobeReplaced = nucleotideSequence.substring(nucleotideSequence.indexOf('('), nucleotideSequence.indexOf(')')+1);
			if(stringTobeReplaced.length()<5){
				nucleotideAfterMutation=nucleotideBeforeMutation=nucleotideSequence.substring(0, nucleotideSequence.indexOf('(')).
						concat(Character.toString(nucleotideSequence.charAt(nucleotideSequence.indexOf('(')+1))).
						concat(nucleotideSequence.substring(nucleotideSequence.indexOf(')')+1));
				System.out.println(nucleotideAfterMutation);
				
			}
			else{
				String mutatedNucleotide = Character.toString(nucleotideSequence.charAt(nucleotideSequence.indexOf('>')+1));
				nucleotideAfterMutation= nucleotideAfterMutation.replace(stringTobeReplaced,mutatedNucleotide);
				String unMutatedNucleotide = Character.toString(nucleotideSequence.charAt(nucleotideSequence.indexOf('>')-1));
				nucleotideBeforeMutation = nucleotideBeforeMutation.replace(stringTobeReplaced,unMutatedNucleotide);
			}
		}
		else if(nucleotideSequence.contains("[")){
			//System.out.println("There has been a mutation of type: Insertion");
			stringTobeReplaced = nucleotideSequence.substring(nucleotideSequence.indexOf('['), nucleotideSequence.indexOf(']')+1);
			nucleotideAfterMutation = (nucleotideBeforeMutation.replace("[","")).replace("]", "");
			nucleotideBeforeMutation = 	nucleotideAfterMutation.replace(stringTobeReplaced, "");	
		}
		else if( nucleotideSequence.contains("(")){
			//System.out.println("There has been a mutation of type: Deletion");
			stringTobeReplaced = nucleotideSequence.substring(nucleotideSequence.indexOf('('), nucleotideSequence.indexOf(')')+1);
			nucleotideBeforeMutation = (nucleotideBeforeMutation.replace("(","")).replace(")", "");	
			nucleotideAfterMutation = 	nucleotideAfterMutation.replace(stringTobeReplaced, "");	
		}
		nucleotideAuthCrpt.put(nucleotideBeforeMutation, nucleotideAfterMutation);
		return nucleotideAuthCrpt;
	}
	public static void writeCrypticSpliceSiteToFile(){
		try{
		    PrintWriter writer = new PrintWriter("CrypticSpliceSite.txt", "UTF-8");
		    for ( String entry : crypticSpliceSites) {
				if(!entry.equalsIgnoreCase("ZZZZZZZZZ")){
					 writer.println(entry);
				}
			}
		    writer.close();
		} catch (IOException e) {
		   e.printStackTrace();
		}
	}
	public static void writeNeighborSpliceSiteToFile(){
		try{
		    PrintWriter writer = new PrintWriter("NeighboringSpliceSite.txt", "UTF-8");
		    for ( String entry : neighboringSpliceSites) {
					 writer.println(entry);
			}
		    writer.close();
		} catch (IOException e) {
		   e.printStackTrace();
		}
	}
	public static void writeWebsitesWithLessNeighbors(String WebsitesWithLessNeighbors){
		BufferedWriter bw = null;
		FileWriter fw = null;
		try {
			File file = new File("WebsitesWithLessNeighbors.txt");
			// if file doesnt exists, then create it
			if (!file.exists()) {
				file.createNewFile();
			}
			// true = append file
			fw = new FileWriter(file.getAbsoluteFile(), true);
			bw = new BufferedWriter(fw);
			bw.write(WebsitesWithLessNeighbors);
			bw.newLine();
		} catch (IOException e) {

			e.printStackTrace();

		}  finally {
			try {
				if (bw != null)
					bw.close();
				if (fw != null)
					fw.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}
	public static void getNeighborData(String afterMutation,String link,String cryptic9Mer,String authentic9Mer){
		String upstream = afterMutation.substring(0, afterMutation.indexOf("/")+1);
		String downstream = afterMutation.substring(afterMutation.indexOf("/"));
		checkAndAddNeighbors(upstream,downstream,link,cryptic9Mer,authentic9Mer);
	}
	public static void checkAndAddNeighbors(String neighbor9merup,String neighbor9merdown,String link,String cryptic9Mer,String authentic9Mer){
		//System.out.println(link);
		String neighbor = "";
		if(neighbor9merup.length()<100){
			neighbor.concat(neighbor9merup);
			writeWebsitesWithLessNeighbors("Upstream less than 100");
			writeWebsitesWithLessNeighbors(link+"+++"+neighbor9merup.length());
		}
		else
		{
			neighbor = (neighbor9merup.substring(neighbor9merup.indexOf("/")-99, neighbor9merup.indexOf("/")));
		}
		if(neighbor9merdown.length()<100){
			neighbor.concat(neighbor9merdown);
			writeWebsitesWithLessNeighbors("DownStream less than 100");
			writeWebsitesWithLessNeighbors(link+"+++"+neighbor9merdown.length());
		}
		else
		{
			neighbor.concat(neighbor9merdown.substring(neighbor9merdown.indexOf("/"), 100));
		}
		for(int i=3;i<neighbor.length()-7;i++){
			if(neighbor.substring(i, i+2).equals("GT")){	
				String neighbor9mer = neighbor.substring(i-3, i+6);
				if(!neighbor9mer.equalsIgnoreCase(cryptic9Mer)&&neighbor9mer.equalsIgnoreCase(authentic9Mer)){
					//System.out.println(link + neighbor9mer +" ++++++++++++++ "+ authentic9Mer+" ++++++++ "+cryptic9Mer);
				}
				else if(!neighbor9mer.equalsIgnoreCase(cryptic9Mer)&&!neighbor9mer.equalsIgnoreCase(authentic9Mer)){
					neighboringSpliceSites.add(neighbor9mer);
				}
			}
		}
	}
	public static void main(String[] args) {
		try{
			WebClient webClient = new WebClient();
			HtmlPage page = webClient.getPage(URL);
			getPageLinks(webClient,page);
			get9mers();
			System.out.println("done");
			writeCrypticSpliceSiteToFile();
			writeNeighborSpliceSiteToFile();
		}catch (FailingHttpStatusCodeException e) {
			
			e.printStackTrace();
		} catch (MalformedURLException e) {
			
			e.printStackTrace();
		} catch (IOException e) {
			
			e.printStackTrace();
		}	
	}
}
