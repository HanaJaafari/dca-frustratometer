#!/usr/bin/python
#Reads Seqs.fa and outputs the alignments
#Writes URLS.txt as a backup to restart jackhmmer if programs breaks
#Everything is hardcoded


def gotoMainPage(browser,Restart_url):
    import time
    while browser.title<>'score results | HMMER' or browser.current_url<>Restart_url:
        print "In %s"%browser.title
        print "Going to main page"
        browser.get(Restart_url)
        time.sleep(5)

def HMMER_onlinev2(name,seq,restart=None,database='Reference Proteomes'):
    '''Does multiple iterations on jackhammer until the same number of sequences
    (or less) is found three succesive for the exact domain. Downloads from the
    set of three iterations that has the maximum number of sequences.'''
    
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
    from selenium.webdriver.support.ui import Select
    import selenium
    import time
    import urllib
    
    #Start the first two iterations
    print "Opening browser"
    browser = webdriver.Firefox()
    browser.implicitly_wait(5)
    #webdriver.manage().timeouts().implicitlyWait()
    if not restart:
        print "Starting first iteration"
        
        print "Going to submission page"
        browser.get('https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer')
        
        
        print "Writing sequence"
        seq_box=browser.find_element_by_xpath(".//*[@id='seq']") #Sequence box
        seq_box.send_keys(seq)
        
        print "Selecting database %s"%database
        select = Select(browser.find_element_by_xpath(".//*[@id='sequencedb']/div[1]/label/select")); #Database selection
        #print [o.text for o in select.options]
        select.select_by_visible_text(database);
        
        
        print "Submitting"
        browser.find_element_by_xpath(".//*[@id='subbutton']").click()
        time.sleep(2) #Wait
        
        print "Starting second iteration"
        for i in range(20):
            print "Pressing button"
            try:
                browser.find_element_by_xpath(".//*[@id='jackhmmer_nav']/p/input").click()
                break #Start second iteration
            except selenium.common.exceptions.NoSuchElementException:
                pass
            time.sleep(1)
        time.sleep(25) #Wait
        
        print "Saving URL"
        with open('URLS.txt','a+') as URLs:
            URLs.write('%s %s %s\n'%(name,browser.current_url,database))
        Restart_url=browser.current_url
        print Restart_url
    
    else:
        print "Restarting previous iteration for %s"%name
        browser.get(restart) #Go to seq URL if saved
        Restart_url=restart

    #Go to page if not there yet    
    gotoMainPage(browser,Restart_url)

    #Look for the exact domain of the query:
    print "Looking for the exact domain"
    browser.find_element_by_xpath(".//*[@id='batch_sum']/tbody/tr[%i]/td[2]/a"%1).click() #Go to First iteration
    browser.find_element_by_xpath(".//*[@id='subnav']/ul/li[3]/a").click() #Go to Domain section
    exact_name=browser.find_element_by_xpath(".//*[@id='exact']/p[2]/strong") #Locate exact domain name
    Original_domain=exact_name.text
    n=browser.find_element_by_xpath(".//*[@id='exact']/a/span[1]")#Locate exact domain number of sequences
    number_sequences=int(n.text.split('\n')[0])
    print "Original domain found: %s"%Original_domain
    print "Number of sequences on original domain: %i" %number_sequences
    
    #Do more iterations until convergence
    Converged=0
    I=1
    for i in range(2,50):
        
        if number_sequences>10000:
            print "The number of sequences is too large, will change the working database from %s"%database
            if database=='Reference Proteomes':
                database='rp75'
            elif database=='rp75':
                database='rp55'
            elif database=='rp55':
                database='rp35'
            elif database=='rp35':
                database='rp15'
            elif database=='rp15':
                print "No more databases avaible"
                return
            print "to %s"%database 
            browser.close()
            HMMER_onlinev2(name,seq,restart=None,database=database)
            return   
        print "On iteration %i"%i
        gotoMainPage(browser,Restart_url)
        print "Checking if iteration ran already"
        try: #Check if there is a new iteration and go
            browser.find_element_by_xpath(".//*[@id='batch_sum']/tbody/tr[%i]/td[2]/a"%i).click() #Item in table
        except selenium.common.exceptions.NoSuchElementException:
            #Try to do another iteration
            print "Iteration not run yet, will submit new iteration"
            for k in range(20):
                #Go to page if not there yet
                gotoMainPage(browser,Restart_url)
                #Do next iteration
                print "Submitting next iteration"
                try:
                    submit_button=browser.find_element_by_xpath(".//*[@id='next_iteration']/input[3]")
                    if submit_button.is_enabled():
                        submit_button.click()
                    time.sleep(2) #Start next iteration
                except selenium.common.exceptions.NoSuchElementException:
                    print "Seems like jackhmmer has converged"
                    Converged=100
                    break
                
                try: #Check if there is a new iteration and go
                    browser.find_element_by_xpath(".//*[@id='batch_sum']/tbody/tr[%i]/td[2]/a"%i).click() #Item in table
                    print "Yeah!, new iteration"
                    break
                except selenium.common.exceptions.NoSuchElementException:
                    print "Iteration has not run yet, trying again"
                    continue
        if Converged>20:
            "No more iterations will be run"
            break
        
        print "Looking for convergence of domain"
        browser.find_element_by_xpath(".//*[@id='subnav']/ul/li[3]/a").click()#Go to Domain section
        #Look if there is an increase on the number of domains
        found=False
        for j in range(1,10):
            print "Reading domain %i"%j
            try:
                domain_name=browser.find_element_by_xpath(".//*[@id='content']/div[5]/ul[2]/li[%i]/p[2]/strong"%j).text
                print "Domain %i: %s"%(j,domain_name)
                if domain_name<>Original_domain:
                    print "Domain %s different from Original domain %s"%(domain_name,Original_domain)
                    continue
                else:
                    print "Domain %s found"%domain_name
                    found=True
                    J=j
                    n=browser.find_element_by_xpath(".//*[@id='content']/div[5]/ul[2]/li[%i]/a/span[1]"%j)
                    n=int(n.text.split('\n')[0])
                    if n<=number_sequences:
                        Converged+=1
                        if n==number_sequences:
                            print "%i sequences found again!"%n
                            I=i
                        else:
                            print "Sequences in domain are now only %i sequences (before %i)"%(n,number_sequences)
                    else:
                        Converged=0
                        print "There are now %i sequences (before %i)"%(n,number_sequences)
                        number_sequences=n
                        I=i
                    break
            except selenium.common.exceptions.NoSuchElementException:
                try:
                    browser.find_element_by_xpath(".//*[@id='wrapper']/div[2]/div[3]/ul/li[%i]/a/span[1]"%j)
                    print "Domain without architecture found"
                    continue
                except selenium.common.exceptions.NoSuchElementException:
                    break
        if not found:
            print "No, Domain dissapeared!"
            if i==2:
                j=1
            else:
                j=J
            i=i-1
            break
        if Converged==1:
            print "Almost Converged..."
        if Converged>1:
            print "Yeah! Converged!"
            break
    #Download
    print "Downloading multi-alignment from iteration %i and domain %s"%(I,Original_domain)
    gotoMainPage(browser,Restart_url)
    browser.find_element_by_xpath(".//*[@id='batch_sum']/tbody/tr[%i]/td[2]/a"%I).click() #Item in table
    browser.find_element_by_xpath(".//*[@id='subnav']/ul/li[3]/a").click()#Domain section
    browser.find_element_by_xpath(".//*[@id='content']/div[5]/ul[2]/li[%i]/p[1]/a"%j).click() #Select correct domain
    time.sleep(10)
    browser.find_element_by_xpath(".//*[@id='subnav']/ul/li[4]/a").click() #Download section
    #browser.find_element_by_xpath(".//*[@id='format']/div[2]/a[9]").click()
    Download_file=browser.current_url+'?format=afa'
    #Download_file=Download_file.replace('/results/','/download/')
    print Download_file
    browser.close()
    urllib.urlretrieve(Download_file, "%s_alignment.fasta.gz"%name)
    

if __name__=='__main__':
    from Bio import SeqIO
    import os
    records = list(SeqIO.parse("Seqs.fa", "fasta"))
    Restart={}
    if os.path.isfile('URLS.txt'): 
        with open('URLS.txt') as URLs:
            for line in URLs:
                a,b=line[:-1].split()[0:2]
                try:
                    c=line[:-1].split()[2]
                except IndexError:
                    c=''
                if c in ['rp75','rp55','rp35','rp15']:
                    Restart.update({a:[b,c]})
                else:
                    Restart.update({a:[b]})
                
    #print Restart
    for i,record in enumerate(records):
        if os.path.isfile("%s_alignment.fasta.gz"%record.name):
            continue
        elif record.name in Restart.keys():
            database='Reference Proteomes'
            if len(Restart[record.name])==2:
                restart=Restart[record.name][0]
                database=Restart[record.name][1]
            else:
                restart=Restart[record.name][0]
            print "Restart page: %s"%restart
            print "Restart database: %s"%database
            HMMER_onlinev2(record.name,str(record.seq),restart,database)
        else:
            HMMER_onlinev2(record.name,str(record.seq))


