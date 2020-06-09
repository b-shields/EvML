# The game app is focused on the CH arylation study from recent BMS Princeton PCI collaboration
# Original game design idea was from Ben Shields from Princeton
# The Shiny App was prepared Jun Li, Jay Stevens, and Jake Janey from BMS and Ben Shields

# Imports

library(shiny)
library(dplyr)
library(ggplot2)
library(reshape2)

# Read data

conditions<-read.csv("All_CH_arylation_Experiments_Jay_12052019.csv", header=TRUE, stringsAsFactors=FALSE,sep = ',')  

conditions[,c(1:5)]<-lapply(conditions[,c(1:5)],as.character)

conditions$yield<-round(conditions$yield,1)

#base     ligand Solvent concentration temperature yield
#KOAc  BrettPhos    DMAc           0.1         105   5.5
#KOAc    PPhtBu2    DMAc           0.1         105   0.0
#KOAc tBPh-CPhos    DMAc           0.1         105  79.0
#KOAc  PCy3 HBF4    DMAc           0.1         105   7.3

all_condition_index<-as.character(seq(1:nrow(conditions))) #used for random selection purpose

#############Assign the percentiles of the yield################
#This is to show the prior participant performances to avoid using the actual optimal yields 
#to inform future players

batches<-5
percentile_table <- conditions %>% 
                mutate(yield_percentile = percent_rank(yield))

#Read all the participants selections from the result folder
filenames <- list.files("results", pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, function(a) read.csv(a,colClasses=c(rep("character",4))))

#Very important: In the stored csv file, the last row is the 1st selection!
#Therefore, to evaluate in leader board, the batch # where the max was observed should be calculated backwards 

userTab<-NULL
cumTab<-NULL
for (i in 1:length(ldf)){
  ldf[[i]]<- ldf[[i]] %>% left_join (percentile_table, by=c('base'='base',
                                                            'ligand'='ligand',
                                                            'solvent'='solvent',
                                                            'concentration'='concentration',
                                                            'temperature'='temperature')) %>%
    mutate(user = strsplit(strsplit(filenames[i],"/")[[1]][2],"_")[[1]][1])
  
  #Prepare a datatable to show where is the highest observed percentile achieved at which 5-run batch by who
  whichRxnBest<-nrow(ldf[[i]])-which.max(ldf[[i]]$yield_percentile)+1  #count backwards from last row since it is the 1st rxn
  best_sofar<-c(ceiling(whichRxnBest/batches),max(ldf[[i]]$yield_percentile),unique(ldf[[i]]$user))
  userTab<-rbind(userTab,best_sofar)
  
  store<-cbind(seq(1:nrow(ldf[[i]])),ldf[[i]][,c("yield_percentile","user")])
  colnames(store)[1]<-"Exp_No"
  cumTab<-rbind(cumTab,store)
}
colnames(userTab)<-c("Batch_No","Observed_Max_Yield_Percentile","User")
rownames(userTab)<-NULL
userTab<-as.data.frame(userTab)
userTab$Observed_Max_Yield_Percentile<-as.numeric(as.character(userTab$Observed_Max_Yield_Percentile))
#userTab_order<-userTab[order(userTab$Observed_Max_Yield_Percentile,rev(userTab$Batch_No),decreasing=T),] #if there is a tie in max yield, the lowest batch will rank on top
userTab$Batch_No<-as.numeric(as.character(userTab$Batch_No))
userTab_order<-dplyr::arrange(userTab,desc(Observed_Max_Yield_Percentile),Batch_No)

cumTab<-as.data.frame(cumTab)
name_list<-toupper(as.character(userTab$User))

###############################################################################
fields <- colnames(conditions)
count<-0
outputDir<-"results"

design_table = data.frame(base = NULL, ligand = NULL,
                          solvent = NULL, concentration = NULL, temperature = NULL)
design_table_history = data.frame(base = NULL, ligand = NULL,
                                  solvent = NULL,  concentration = NULL, temperature = NULL)
design_table_current = data.frame(Base = NULL, Ligand = NULL,
                                  solvent = NULL,  concentration = NULL, temperature = NULL)
retrieve_table = data.frame(base = NULL, ligand = NULL,
                            solvent = NULL,  concentration = NULL, temperature = NULL, yield=NULL)
random_table_history = data.frame(base = NULL, ligand = NULL,
                                  solvent = NULL,  concentration = NULL, temperature = NULL,yield=NULL) #need incl yield in random table


ui<-shinyUI(
  navbarPage(
    "EvML: Reaction Optimization",
    tabPanel("Information",
             fluidPage( withMathJax(includeMarkdown("README.md")))
    ),
    
    tabPanel("Run Experiments",
             sidebarPanel(
               helpText("Please enter your user name first to initiate the game."),
               textInput("nname", "Name", value = "", width = NULL, placeholder="Enter name here to begin"),
               selectInput("human","ID",c("Faculty","Graduate Student","Postdoc","Med Chemist","Process Chemist","Engineer","Statistician","Other")),
               selectInput("type","Type",c("Academic", "Pharma","Biotech","CRO or CMO","Other")),
               selectInput("year","Experiences in Pd cross-coupling",c("1~5 years","5~10 years","Above 10 years","No experience")),
               br(),
               helpText("Reactions:"),
               selectInput("base","Base",unique(conditions$base)),
               selectInput("ligand","Ligand",unique(conditions$ligand)),
               selectInput("solvent","Solvent",unique(conditions$solvent)),
               selectInput("concentration","Concentration",unique(conditions$concentration)),
               selectInput("temperature","Temperature",unique(conditions$temperature)),
               conditionalPanel( condition = "output.modl != 0 || output.count == 0 || output.status == 'Run'",
                                 actionButton("add","Add")),

               actionButton("delete_btn", "Delete"),
               helpText("Select 5 conditions to run in parallel per batch. If you need to delete a condition prior to submission, highlight the row and hit delete. The 'Run experiments' button will show up when 5 experiments have been selected. Note: please try to avoid selecting duplicate conditions; if you do an error message will be displayed."),
               br(),
               
               #make sure the max output count no larger than 100 runs below, which are five-run per batch and 20 batches
               conditionalPanel( condition = "output.modl == 0 && output.count != 0 && output.count <= 100 && output.status != 'Run'",
                                 actionButton("retrieve","Run experiments"))
                               
             ),
             
             # a table of reactive outputs
             mainPanel(
               mainPanel(
                 textOutput("message"),
                 tags$head(tags$style("#message{color: red;
                                 font-size: 40px;
                                      font-style: bold;
                                      }"
                         )
                 ),
                 img(src="img.png", height="100%", width="120%", align = "middle"),
                 DT::dataTableOutput("responses", width = 500), 
                 tags$hr(),
                 DT::dataTableOutput("responsesCurrent",width = 500)
               )
             )
    ),
    
    tabPanel("Experiment Results",
             
             mainPanel(
               mainPanel(
              #   img(src='table.PNG', height="100%", width="150%", align = "middle"),
                 DT::dataTableOutput("responses2", width = 1000), tags$hr(),
                 plotOutput("plot2")
               )
             )
    ),
    
    tabPanel("Leader Board",
             
             mainPanel(
               mainPanel(
                 DT::dataTableOutput("leader_board", width = 1000), tags$hr(),
                 helpText("Above are the results from other users by percentile. There are no live results here."),
                 plotOutput("plot3")
               )
             )
    )
  )
)
  
server<-shinyServer(function(input, output, session){  
    
    counter2 <- reactiveValues(status="NoRun")
    counter <- reactiveValues(countervalue = 0)
    values <- reactiveValues(fileName = NULL)
    
    design_table_history <- reactiveVal(design_table_history) #use to lock the previous rows to distable delete
    design_table_current <- reactiveVal(design_table_current)
    random_table_history <- reactiveVal(random_table_history)
    retrieve_table <- reactiveVal(retrieve_table)

    loadData2 <- reactive({
      #Caveats: if there are duplicates in the design dataset
      #the merge function whether using sort=FALSE or not will cause the display order 
      #change in yield table. The original mapply work w/o issue but ugly.
      #The dplyr left_join has an issue to show the 1st empty set in yield section!
      #Due to the fact thet the actual order incl duplicates will be stored in csv for later analysis,
      #the display in yield table is just for user check up only.  Therefore, we kept 
      #the merge function here and also apply unique view to filter the duplicates
      
      unique(merge(retrieve_table(),conditions,sort=FALSE)) 
    })

    Warn<-reactive({
      validate(
        need(
          !toupper(input$nname) %in% name_list, 
          'Same nickname used before, Try new one'
        ),
      ''
      )
    })
    
    getID <- reactive({
      id <- sapply("human", function(x) input[[x]])
      id
    })

    getYear <- reactive({
      experience <- sapply("year", function(x) input[[x]])
      experience
    })
    
    getType <- reactive({
      ind <- sapply("type", function(x) input[[x]])
      ind
    })
    
    #Make sure the new nickname does not exist before
    getNickname <- reactive({
      id <- sapply("nname", function(x) input[[x]])
      id
    })

    observeEvent(input$add, {
        counter2$status <- "NoRun"
        counter$countervalue <- counter$countervalue + 1  #For redundant choice and unresponsive, extra clicks were recorded!
        current <-data.frame(base = input$base,
                             ligand = input$ligand, 
                             solvent = input$solvent,
                             concentration = input$concentration,
                             temperature = input$temperature) #each added current condition
        total_current <- rbind(current,design_table_current()) #whole block of added current conditions
        design_table_current(total_current)   #store in a reactive container for future delete
 
        nickname2<-ifelse(getNickname()=="","NoName",getNickname())  
        values$fileName<- sprintf("%s_%s_%s_%s.csv", nickname2,getType(),getID(),getYear())
    })

    observeEvent(input$retrieve, {
      counter2$status <- "Run"
      tmp <- design_table_current() #store the modifed most current block
      tmp2 <- rbind(tmp,design_table_history()) #combine with previous history runs
      design_table_history(tmp2)
      
      write.table(x = design_table_history(),
                  file.path(outputDir, values$fileName),
                  append=FALSE,col.names=TRUE,sep = ",",row.names=FALSE)  #save the updated history
      
      #after click this run/retrieve section and update the csv file
      #we need to clean the reactive design_table_current memory
      #therefore, after run, there are no rows to delete!
      update_current = data.frame(base = NULL, ligand = NULL,
                                        solvent = NULL, concentration = NULL, temperature = NULL)
      design_table_current(update_current)
      retrieve_table(tmp2)  #All historical data to get the yield and output
      
      #random table/ just to match the actual layout where new conditions selected on top of the old ones
      remain_index <- setdiff(all_condition_index,rownames(random_table_history()))
      tmp_rnd <- rbind(conditions[as.numeric(sample(remain_index,batches,replace=FALSE)),],
                       random_table_history())
      random_table_history(tmp_rnd)
    })
    
    observeEvent(input$delete_btn, {  #This only allows to delete most current working table
      counter$countervalue <- counter$countervalue - 1
      total = design_table_current()
      #print(nrow(total))
      if (!is.null(input$responsesCurrent_rows_selected)) {
        total <- total[-as.numeric(input$responsesCurrent_rows_selected),]
      }
      design_table_current(total)
    })
    
    yield_extracted <- reactive({
      loadData2() %>%
                select (yield) %>%
                mutate(id = row_number()) %>%
                     cbind(yield_random_robot=random_table_history()$yield) %>%
                        melt(id.vars = "id")
    })
    
    
    # Show the previous responses
    # (update with current response when add is clicked)
    output$responses <- DT::renderDataTable({
      DT::datatable(
        design_table_history(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          lengthChange = FALSE
        )
    )})

    output$responsesCurrent <- DT::renderDataTable({
      DT::datatable(
        design_table_current(),
        selection = 'single',
        rownames = FALSE,
        options = list(
          #dom = 't',
          pageLength = 5,
          lengthChange = FALSE
        )
      )})
    
    #yield tab output
    output$responses2 <- DT::renderDataTable({input$retrieve   
      DT::datatable(
        loadData2(),  #reactive data
        rownames = FALSE,
        options = list(
          pageLength = 10,
          lengthChange = FALSE
        )
      )})

    output$message <- renderText({
      Warn()
    })
    
    
    #Leader board output, not dynamic, just show the historical ones
    output$leader_board <- DT::renderDataTable(userTab_order,   
        options = list(
          pageLength = 10,
          lengthChange = FALSE
        ),
        rownames= FALSE
      ) 
        
    output$plot2<-renderPlot({
      input$retrieve
      ggplot(yield_extracted(),aes(x=id,y=value,color=variable)) +
             geom_point(size=6) +
             geom_line() +
             xlab("Sequence number") + ylab("Yield") +
             ggtitle("Progress on Reaction Condition Exploration") + 
             theme(plot.title = element_text(size=16, face="bold", vjust=2))},
      height = 400,width = 1000)
    
    output$plot3<-renderPlot({
      ggplot(cumTab,aes(x=Exp_No,y=yield_percentile,color=user)) +
        geom_point(size=6) +
        geom_line() +
        xlab("Experiments") + ylab("Observed Yield by Percentile") +
        ggtitle("Comparison of Reaction Performance by Users") + 
        theme(plot.title = element_text(size=16, face="bold", vjust=2))
      },
      height = 400,width = 1000)
    
    output$count <- reactive({
          counter$countervalue
    })
    
    output$modl <- reactive({
      counter$countervalue %% batches       #Here is 5 runs per batch for screen lock purpose
    })
    
    output$status <- reactive({
      counter2$status
    })
    
    outputOptions(output, "modl", suspendWhenHidden = FALSE)
    outputOptions(output, "count", suspendWhenHidden = FALSE)
    outputOptions(output, "status", suspendWhenHidden = FALSE)
})


shinyApp(ui, server)  
