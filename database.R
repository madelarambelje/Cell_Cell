# This script will create a custom DataBase for CellChat, based on ICELLNET's DB
# Nick Veltmaat
# 16-5-2021

if (DB == "CellChat") {
  # Set CellChat database / reference
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  cellchat1@DB <- CellChatDB    #Set DB for Healthy Subset 1
  cellchat2@DB <- CellChatDB   #Set DB for Inflamed Subset 2
  
} else if (DB  == "ICELLNET") {
  # First Load the ICELLNETDB, which is stored as an RDS file.
  ICELLNETDB <- readRDS('Vignette/data/Icellnet_db.rds')
  
  CellChatDB <- CellChatDB.human # Load CellChatDB
  # Subset DB based on ICELLNET's DB
  CellChatDB$interaction <- CellChatDB$interaction[CellChatDB$interaction$ligand %in% ICELLNETDB$`Ligand 1`, ]
  CellChatDB$complex <- CellChatDB$complex[CellChatDB$complex$subunit_2 %in% ICELLNETDB$`Receptor 2`, ]
  showDatabaseCategory(CellChatDB)
  cellchat1@DB <- CellChatDB    #Set DB for Healthy Subset 1
  cellchat2@DB <- CellChatDB   #Set DB for Inflamed Subset 2
  
} else {
  print("Database not found, try setting DB to 'ICELLNET' or 'CellChat'")
}
