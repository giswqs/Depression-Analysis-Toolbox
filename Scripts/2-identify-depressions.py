__author__ = 'Qiusheng'

import arcpy
import numpy
import os
from arcpy import env
import datetime

def PreProcess(DEMRasterPath,min_size,buffer_dist):

    arcpy.CheckOutExtension("Spatial")
    input_dir = os.path.split(DEMRasterPath)[0]
    input_ras = os.path.split(DEMRasterPath)[1]
    #min_size = 100
    #buffer_dist = 5

    env.workspace = input_dir
    env.overwriteOutput = True

    if arcpy.Exists(input_ras) == False:
        print("The input raster does not exist")
        quit()


    ### Mean Focal Statistics
    ras_mf = arcpy.sa.FocalStatistics(input_ras,"Rectangle 3 3 CELL","MEAN","DATA")
    #mean_filter_file = os.path.join(input_dir,"ras_mf")
    #ras_mf.save("ras_mf")

    ### Fill depression
    ras_fill = arcpy.sa.Fill(ras_mf)
    #ras_fill.save("ras_fill")

    ### Get sink
    ras_sink = ras_fill - ras_mf
    #ras_sink.save("ras_sink")

    ### Convert foreground sink
    #ras_sink_fg = arcpy.sa.Con(ras_sink,ras_sink,"#",'"Value">0')
    #ras_sink_fg.save("ras_sink_fg")

    ### Convert sink to binary image
    ras_sink_bin = arcpy.sa.Con(ras_sink > 0,1)
    #ras_sink_bin.save("ras_sink_bin")

    ### Region group
    ras_region = arcpy.sa.RegionGroup(ras_sink_bin,"FOUR","WITHIN","ADD_LINK")
   # ras_region.save("ras_region")

    ### Convert raster to polygon
    region_poly_name = os.path.join(input_dir,"region_poly.shp")
    arcpy.RasterToPolygon_conversion(ras_region,region_poly_name,"NO_SIMPLIFY")

    ### Select polygon based on minimum size
    ##input_dir = r"D:\Dropbox\Research\Codes\Python\SurfaceDepression\data\testdata.gdb"
    area_field = "Area"
    arcpy.AddField_management(region_poly_name, area_field, "DOUBLE")
    arcpy.CalculateField_management(region_poly_name,"Area","!shape.area@squaremeters!","PYTHON_9.3","#")
    sqlExp = area_field + ">=" + str(min_size)

    ##print sqlExp
    region_poly_select_name = os.path.join(input_dir,"region_poly_select.shp")
    arcpy.Select_analysis(region_poly_name,region_poly_select_name,sqlExp)

    ### Calculate field
    arcpy.CalculateField_management(region_poly_select_name,"gridcode","1","PYTHON")

    #### Convert polygon to raster
    region_poly_ras = os.path.join(input_dir,"region_poly_ras.tif")
    arcpy.PolygonToRaster_conversion(region_poly_select_name,"gridcode",region_poly_ras,"CELL_CENTER","NONE","1")

    ### Convert foreground sink TO 0
    ras_sink_bg = ras_mf - ras_mf
    ras_sink_bg.save(os.path.join(input_dir,"ras_sink_bg"))

    #ras_sink_final = "ras_sink_final"
    in_ras_list = [region_poly_ras,ras_sink_bg]
    ras_sink_final_name = arcpy.sa.CellStatistics(in_ras_list,"SUM","DATA")
    #ras_sink_final_name.save(ras_sink_final)

    ### Convert foreground sink
    dem_name = arcpy.sa.Con(ras_sink_final_name==1,ras_mf,ras_fill)
    dem_name.save(os.path.join(input_dir,"dem_final.tif"))

    ### buffer sink area
    sink_buffer_name = os.path.join(input_dir,"sink_buffer.shp")
    sqlExp = str(buffer_dist) + " Meters"
    arcpy.Buffer_analysis(region_poly_select_name,sink_buffer_name,sqlExp,"","","ALL","")

    dem_sink = arcpy.sa.ExtractByMask(dem_name,sink_buffer_name)
    dem_sink.save(os.path.join(input_dir,"dem_sink_buffer.tif"))

    dem_sink_depth = ras_fill - dem_name
    dem_sink_depth_name = arcpy.sa.Con(dem_sink_depth>0,dem_sink_depth)
    dem_sink_depth_name.save(os.path.join(input_dir,"dem_sink_depth.tif"))



    arcpy.Delete_management(region_poly_name)
    arcpy.Delete_management(region_poly_select_name)
    arcpy.Delete_management(region_poly_ras)
    arcpy.Delete_management(sink_buffer_name)
    arcpy.Delete_management(os.path.join(input_dir,"dem_sink_buffer.tif"))
    arcpy.Delete_management(os.path.join(input_dir,"ras_sink_bg"))
    arcpy.Delete_management(os.path.join(input_dir,"dem_final.tif"))


    return os.path.join(input_dir,"dem_sink_depth.tif")

    print "PreProcess Done!"



def JoinAllTables(inContour_Poly_ID,inStastic,inSummarize, inContour_Polyline):
    arcpy.AddIndex_management (inContour_Polyline, "Cont_ID", "Cont_Index")
    arcpy.AddIndex_management (inContour_Poly_ID, "Cont_ID", "Cont_Index")
    arcpy.AddIndex_management (inStastic, "Cont_ID", "Cont_Index")
    if (os.path.splitext(inSummarize)[1].lower() == ".dbf"):
        arcpy.AddIndex_management (inSummarize, "src_Cont_I", "Cont_Index")
        arcpy.JoinField_management(inContour_Poly_ID,"Cont_ID",inSummarize,"src_Cont_I")
    else:
        arcpy.AddIndex_management (inSummarize, "src_Cont_ID", "Cont_Index")
        arcpy.JoinField_management(inContour_Poly_ID,"Cont_ID",inSummarize,"src_Cont_ID")

    arcpy.JoinField_management(inContour_Poly_ID,"Cont_ID",inStastic,"Cont_ID")
    arcpy.JoinField_management(inContour_Poly_ID,"Out_Nbr_ID",inContour_Polyline,"Cont_ID")
    #arcpy.AddJoin_management(inContour_Poly_ID,"Cont_ID",inPolygon_Nei,"Cont_ID","KEEP_ALL")
    arcpy.AddMessage("JoinAllTables Over!")

def UpdateAllNbrInfo(inContour_Poly_ID,inPoly_Nei):
    AllInfoCount = len(inContour_Poly_ID)
    for polygon_ID_Index in range(0,AllInfoCount):
        self_ID = inContour_Poly_ID[polygon_ID_Index][0]

        if (os.path.splitext(inPoly_Nei)[1].lower() == ".dbf"):
            selection = "src_Cont_I = " + str(self_ID)
            cursor = arcpy.SearchCursor(inPoly_Nei,selection,fields="src_Area;nbr_Area;nbr_Cont_I",sort_fields="nbr_Area D")
        else:
            selection = "src_Cont_ID = " + str(self_ID)
            cursor = arcpy.SearchCursor(inPoly_Nei,selection,fields="src_Area;nbr_Area;nbr_Cont_ID",sort_fields="nbr_Area D")
        record = cursor.next()
        try:
            if record:
                self_area = float(record.getValue("src_Area"))
                nbr_area = float(record.getValue("nbr_Area"))
                if(nbr_area > self_area):
                    if (os.path.splitext(inPoly_Nei)[1].lower() == ".dbf"):
                        nbr_ID = int(record.getValue("nbr_Cont_I"))
                    else:
                        nbr_ID = int(record.getValue("nbr_Cont_ID"))
                    nbr_poly_tup = ReadPolyInfo(inContour_Poly_ID,nbr_ID)
                    UpdateOutNerInfo(inContour_Poly_ID,self_ID,nbr_ID,int(nbr_poly_tup[6]))
        except:
            a = 1
        del cursor
    arcpy.AddMessage("UpdateAllNbrInfo Over!")

def UpdateOutNerInfo(inContour_Poly_ID,Poly_ID,Out_Nbr_ID,Out_Nbr_Cnt):
    self_ID = str(Poly_ID)
    self_poly = ReadPolyInfo(inContour_Poly_ID,self_ID)
    self_poly[7] = Out_Nbr_ID
    self_poly[8] = Out_Nbr_Cnt


def ReadAllInfo(inContour_Poly_ID):
    if (os.path.splitext(inContour_Poly_ID)[1].lower() == ".shp"):
        Allinfo = arcpy.da.TableToNumPyArray(inContour_Poly_ID,('Cont_ID','FIRST_Cont','SUM_Lap_Ar','MEAN','Line_Cou','Level_ID','NBR_CNT','Out_Nbr_ID','Out_Nbr_Ct','MIN','MAX','LabelFlag','Level_BD'))
    else:
        Allinfo = arcpy.da.TableToNumPyArray(inContour_Poly_ID,('Cont_ID','FIRST_Cont_','SUM_Lap_Area','MEAN','Line_Cou','Level_ID','NBR_CNT','Out_Nbr_ID','Out_Nbr_Ct','MIN','MAX','LabelFlag','Level_BD'),"#",skip_nulls=False,null_value = 0)
    return Allinfo

#0:src_ID
#1:src_elev
#2:src_area
#3:src_mean
#4:src_line_count
#5:src_level
#6:src_nbr_cnt
#7:out_nbr_ID
#8:out_nbr_cnt
#9:src_min_elev
#10:src_max_elev
#11:src_label_flag
#12:src_level_bound
def ReadPolyInfo(Allinfo,Poly_ID):
    src_ID = int(Poly_ID)
    polygon_Index = numpy.where(Allinfo['Cont_ID'] == src_ID)
    polygon = Allinfo[polygon_Index[0][0]]
    return polygon

def CheckOutNeiPoly(Allinfo,self_ID):
    self_poly = ReadPolyInfo(Allinfo,self_ID)
    if int(self_poly[8]) > 0:
        nbr_poly_tup = ReadPolyInfo(Allinfo,int(self_poly[7]))
        return nbr_poly_tup
    else:
        return 0

def LabelContLevel(Allinfo,Poly_ID,level_id):
    self_ID = str(Poly_ID)
    if self_ID == "9847":
        print self_ID
    self_poly = ReadPolyInfo(Allinfo,self_ID)
    self_poly[5] = int(level_id)
    self_poly[11] = int(1)

def LabelBoundLevel(Allinfo,Poly_ID,level_id):
    self_ID = str(Poly_ID)
    if self_ID == "9847":
        print self_ID
    self_poly = ReadPolyInfo(Allinfo,self_ID)
    self_poly[12] = int(level_id)
    self_poly[11] = int(1)

def MinusNbrCount(Allinfo,Poly_ID):
    self_ID = str(Poly_ID)
    self_poly = ReadPolyInfo(Allinfo,self_ID)
    self_poly[6] =  self_poly[6] - 1


PeakAL = []
def IdentifyLevel0(Allinfo):
    # polygons = arcpy.UpdateCursor(inContour_Poly_ID)
    # for polygon in polygons:
    #     ToLevelID.append(polygon.getValue("Cont_ID"))
    todoList = numpy.where((Allinfo['Line_Cou']==1) & (Allinfo['Level_ID'] == 0))
    for polygon_ID_Index in range(0,len(todoList[0])):
        src_polygon = Allinfo[todoList[0][polygon_ID_Index]]
        src_ID = int(src_polygon[0])
        #src_polygon = ReadPolyInfo(Allinfo,src_ID)
        mean_elev = src_polygon[3]
        outline_elev = src_polygon[1]
        if mean_elev <= outline_elev:
            LabelContLevel(Allinfo,src_ID,1)
        else:
            PeakAL.append(src_ID)
            iSaddleId = FindSaddle(Allinfo,src_ID)
            if not(iSaddleId == 0):
                saddle_poly = ReadPolyInfo(Allinfo,iSaddleId)
                if (int(saddle_poly[5]) == 0):
                    LabelContLevel(Allinfo,iSaddleId,1)
    arcpy.AddMessage("IdentifyLevel0 Over!")

def IdentifyLevel1(Allinfo,gMinDepth,gMinArea):
    level_id = 1
    LevelTmpAL = []

    todoList = numpy.where((Allinfo['Level_ID']==1))

    for polygon_ID_Index in range(0,len(todoList[0])):
        src_polygon = Allinfo[todoList[0][polygon_ID_Index]]
        src_ID = int(src_polygon[0])
        iUpperBoundID = FindSameLevelUpperBound(Allinfo,src_ID,level_id)
        if iUpperBoundID > 0:
            temp_poly = ReadPolyInfo(Allinfo,iUpperBoundID)
            max_depth = float(temp_poly[10]) - float(temp_poly[9])
            if max_depth > gMinDepth and float(temp_poly[2]) > gMinArea:
                LabelContLevel(Allinfo,iUpperBoundID,level_id)
                if not (iUpperBoundID in LevelTmpAL):
                    LevelTmpAL.append(iUpperBoundID)
            else:
                if int(temp_poly[7]) > 0:
                    MinusNbrCount(Allinfo,int(temp_poly[7]))
                    LevelTmpAL.append(int(temp_poly[7]))

    for levelTmp in LevelTmpAL:
        CPolyLevel1Tmp = ReadPolyInfo(Allinfo,levelTmp)
        iNbrCnt = int(CPolyLevel1Tmp[6])
        dNbrOutlineElev = float(CPolyLevel1Tmp[1])
        dNbrMeanElev = float(CPolyLevel1Tmp[3])
        dNbr_out_nbr_ID = int(CPolyLevel1Tmp[7])
        if dNbrOutlineElev >= dNbrMeanElev:
            if iNbrCnt == 0:
                LabelContLevel(Allinfo,levelTmp,level_id)
                LabelBoundLevel(Allinfo,levelTmp,level_id)
            elif iNbrCnt == 1 and dNbr_out_nbr_ID == 0:
                LabelContLevel(Allinfo,levelTmp,level_id)
                LabelBoundLevel(Allinfo,levelTmp,level_id)
            elif iNbrCnt == 1 and dNbr_out_nbr_ID > 0:
                #out_nbr_ID = int(CPolyLevel1Tmp[7])
                CPolyLevel1Tmp_Nbr = ReadPolyInfo(Allinfo,dNbr_out_nbr_ID)
                iUpperPolyId = FindSameLevelUpperBound(Allinfo,levelTmp,level_id)

                #CUpperPoly = ReadPolyInfo(inContour_Poly_ID,int(iUpperPolyId))
                LabelContLevel(Allinfo,iUpperPolyId,level_id)
                if int(CPolyLevel1Tmp_Nbr[5]) < level_id:
                    LabelBoundLevel(Allinfo,iUpperPolyId,level_id)
            elif iNbrCnt == 2 and  dNbr_out_nbr_ID > 0:
                iUpperPolyId = FindSameLevelUpperBound(Allinfo,levelTmp,level_id)
                LabelContLevel(Allinfo,iUpperPolyId,level_id)
                LabelBoundLevel(Allinfo,iUpperPolyId,level_id)
    arcpy.AddMessage("IdentifyLevel1 Over!")

def IdentifyLevel2(inContour_Poly_ID):
    level_id = 2
    todoList = numpy.where((inContour_Poly_ID['Level_BD']==1))

    for polygon_ID_Index in range(0,len(todoList[0])):
        src_polygon = inContour_Poly_ID[todoList[0][polygon_ID_Index]]
        src_ID = int(src_polygon[0])
        iUpperLevelPolyId = FindUpperLevel(inContour_Poly_ID,src_ID,level_id)
        if iUpperLevelPolyId > 0:
            LabelContLevel(inContour_Poly_ID,iUpperLevelPolyId,level_id)
    arcpy.AddMessage("IdentifyLevel2 Over!")
def IdentifyAllLevel(inContour_Poly_ID):
    level_id = 2
    while (True):
        todoList = numpy.where((inContour_Poly_ID['Level_ID']==level_id))
        if len(todoList) == 1 and len(todoList[0]) == 0:
            break
        else:
            UpperBoundAL = []
            UpperLevelAL = []
            for polygon_ID_Index in range(0,len(todoList[0])):
                src_polygon = inContour_Poly_ID[todoList[0][polygon_ID_Index]]
                src_ID = int(src_polygon[0])
                iUpperBoundID = FindSameLevelUpperBound(inContour_Poly_ID,src_ID,level_id)
                if iUpperBoundID>0:
                    LabelBoundLevel(inContour_Poly_ID,iUpperBoundID,level_id)
                    UpperBoundAL.append(iUpperBoundID)
            if len(UpperBoundAL) > 0:
                todoList2 = numpy.where((inContour_Poly_ID['Level_BD']==level_id))
                level_id = level_id + 1
                for polygon_ID_Index in range(0,len(todoList2[0])):
                    src_polygon = inContour_Poly_ID[todoList2[0][polygon_ID_Index]]
                    src_ID = int(src_polygon[0])
                    iUpperLevelPolyId = FindUpperLevel(inContour_Poly_ID,src_ID,level_id)
                    if iUpperLevelPolyId > 0:
                        LabelContLevel(inContour_Poly_ID,iUpperLevelPolyId,level_id)
                        UpperLevelAL.append(iUpperBoundID)
                if len(UpperLevelAL) == 0:
                    break
    return level_id
    arcpy.AddMessage("IdentifyLevelAll Over!")

def IdentifyPeakLevel(inContour_Poly_ID):
    for peak in PeakAL:
        peak_poly = ReadPolyInfo(inContour_Poly_ID,int(peak))
        peakNbrAl = []
        peakNbrAl.append(int(peak))
        if int(peak_poly[7]) > 0:
            peak_nbr_poly = ReadPolyInfo(inContour_Poly_ID,int(peak_poly[7]))
            if int(peak_nbr_poly[5]) < 1:
                peakNbrAl.append(int(peak_nbr_poly[0]))
                while(True):
                    peak_poly = peak_nbr_poly
                    if (int(peak_poly[7]) > 0):
                        peak_nbr_poly = ReadPolyInfo(inContour_Poly_ID,int(peak_poly[7]))
                        if int(peak_nbr_poly[5]) < 1:
                            peakNbrAl.append(int(peak_nbr_poly[0]))
                        else:
                            break
                    else:
                        break

        for peakIndex in range(0,len(peakNbrAl)):
            peak_ID = peakNbrAl[int(peakIndex)]
            c_peak_poly = ReadPolyInfo(inContour_Poly_ID,peak_ID)
            peak_nbr_ID = int(c_peak_poly[7])
            if peak_nbr_ID > 0:
                c_peak_nbr_poly = ReadPolyInfo(inContour_Poly_ID,peak_nbr_ID)
                if int(c_peak_nbr_poly[5]) > 0 and float(c_peak_nbr_poly[1]) > float(c_peak_poly[1]):
                    LabelContLevel(inContour_Poly_ID,peak_ID,int(c_peak_nbr_poly[5]))
                else:
                    while (int(c_peak_nbr_poly[7])>0):
                        peak_nbr_ID = int(c_peak_nbr_poly[7])
                        c_peak_nbr_poly = ReadPolyInfo(inContour_Poly_ID,peak_nbr_ID)
                        if int(c_peak_nbr_poly[5]) > 0 and float(c_peak_nbr_poly[1]) > float(c_peak_poly[1]):
                            LabelContLevel(inContour_Poly_ID,peak_ID,int(c_peak_nbr_poly[5]))
                            break
                        else:
                            if int(c_peak_nbr_poly[7])>0:
                                c_peak_nbr_poly = ReadPolyInfo(inContour_Poly_ID,int(c_peak_nbr_poly[7]))
    arcpy.AddMessage("IdentifyLevelPeak Over!")

def PostProcessing(inContour_Poly_ID,gMinDepth,gMinArea):
    AllInfoCount = len(inContour_Poly_ID)
    for polygon_ID_Index in range(0,AllInfoCount):
        self_ID = inContour_Poly_ID[polygon_ID_Index][0]
        self_poly = ReadPolyInfo(inContour_Poly_ID,self_ID)
        self_level = int(self_poly[5])

        if self_level == 0:
            self_out_nbr_ID = int(self_poly[7])
            if self_out_nbr_ID > 0:
                nbr_poly = ReadPolyInfo(inContour_Poly_ID,self_out_nbr_ID)
                if int(nbr_poly[5]) > 0 and float(nbr_poly[1]) >= float(self_poly[1]):
                    LabelContLevel(inContour_Poly_ID,self_ID,int(nbr_poly[5]))
                LabelBoundLevel(inContour_Poly_ID,self_ID,0)
        elif self_level > 0:
            self_out_nbr_ID = int(self_poly[7])
            if self_out_nbr_ID > 0:
                nbr_poly = ReadPolyInfo(inContour_Poly_ID,self_out_nbr_ID)
                if int(nbr_poly[5]) >  int(self_poly[5]):
                    LabelBoundLevel(inContour_Poly_ID,self_ID,self_level)
                elif int(nbr_poly[5]) == int(self_poly[5]):
                    LabelBoundLevel(inContour_Poly_ID,self_ID,0)
                elif int(nbr_poly[5]) == 0:
                    LabelBoundLevel(inContour_Poly_ID,self_ID,self_level)
                elif int(nbr_poly[5]) < int(self_poly[5]):
                    LabelBoundLevel(inContour_Poly_ID,self_out_nbr_ID,self_level)
            else:
                max_depth = float(self_poly[10]) - float(self_poly[9])
                if max_depth <=0:
                    max_depth = float(-1) * max_depth
                if max_depth > gMinDepth and float(self_poly[2]) > gMinArea:
                    LabelBoundLevel(inContour_Poly_ID,self_ID,self_level)
                else:
                    LabelBoundLevel(inContour_Poly_ID,self_ID,0)
    arcpy.AddMessage("postProcess Over!")

def CalculateVolume(inContour_Poly_ID,inDEMRaster):

    raster = arcpy.Raster(inDEMRaster)
    cellSizeHeight = raster.meanCellHeight
    cellSizeWidth = raster.meanCellWidth
    cellArea = cellSizeHeight * cellSizeWidth
    del raster

    if (os.path.splitext(inContour_Poly_ID)[1].lower() == ".shp"):
        outStatisticalTable = os.path.join(os.path.split(inContour_Poly_ID)[0],"tempStastic.dbf")
        outZoneRaster = os.path.join(os.path.split(inContour_Poly_ID)[0],"zoneRaster.tif")
    else:
        database_space = os.path.split(os.path.split(inContour_Poly_ID)[0])[0]
        outStatisticalTable = os.path.join(database_space,"tempStastic")
        outZoneRaster = os.path.join(database_space,"zoneRaster")
    arcpy.CheckOutExtension("Spatial")
    #arcpy.sa.ZonalStatisticsAsTable(inContour_Poly_ID,"Cont_ID",inDEMRaster,outStatisticalTable,"DATA","ALL")
    if os.path.exists(outZoneRaster):
        arcpy.Delete_management(outZoneRaster)
    arcpy.PolygonToRaster_conversion(inContour_Poly_ID,"Cont_ID",outZoneRaster,"CELL_CENTER","",cellSizeHeight)
    arcpy.sa.ZonalStatisticsAsTable(outZoneRaster,"Value",inDEMRaster,outStatisticalTable,"DATA","ALL")

    arcpy.AddField_management(outStatisticalTable, "Cont_ID", "LONG")
    arcpy.CalculateField_management(outStatisticalTable,"Cont_ID","!Value!","PYTHON_9.3")


    arcpy.JoinField_management(inContour_Poly_ID,"Cont_ID",outStatisticalTable,"Cont_ID")
    arcpy.AddField_management(inContour_Poly_ID,"Volume","Double")
    express = "MySub(!COUNT!,!SUM!,!Cont_!) *" + str(cellArea)
    codeblock = """def MySub(count,sum,elev):
        vol = float(elev) * float(count) - float(sum)
        return vol
    """
    arcpy.CalculateField_management(inContour_Poly_ID,"Volume",express,"PYTHON_9.3",codeblock)
    arcpy.Delete_management(outStatisticalTable)
    arcpy.Delete_management(outZoneRaster)





def FindSaddle(inContour_Poly_ID,Poly_ID):
    while(True):
        src_poly = ReadPolyInfo(inContour_Poly_ID,Poly_ID)
        nbr_poly = CheckOutNeiPoly(inContour_Poly_ID,Poly_ID)

        if nbr_poly == 0:
            return 0;
        else:
            nbr_nbr_poly = CheckOutNeiPoly(inContour_Poly_ID,int(nbr_poly[0]))
            if nbr_nbr_poly == 0:
                if float(nbr_poly[1]) > float(src_poly[1]):
                    return int(nbr_poly[0])
                else:
                    return 0
            else:
                nbr_nbr_nbr_poly = CheckOutNeiPoly(inContour_Poly_ID,int(nbr_nbr_poly[0]))
                if nbr_nbr_nbr_poly == 0:
                     if float(nbr_nbr_poly[1]) > float(nbr_poly[1]):
                        return int(nbr_poly[0])
                     else:
                        return 0
                else:
                    if float(src_poly[1]) >=  float(nbr_poly[1]):
                        Poly_ID = int(nbr_poly[0])
                    else:
                        return int(nbr_poly[0])

def FindSameLevelUpperBound(inContour_Poly_ID,Poly_ID,level_id):
    nbr_poly = CheckOutNeiPoly(inContour_Poly_ID,Poly_ID)
    self_poly = ReadPolyInfo(inContour_Poly_ID,Poly_ID)
    out_nbr_cnt = int(self_poly[8])
    if nbr_poly == 0:
        LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
        return Poly_ID

    if (float(nbr_poly[1]) <= float(self_poly[1])):
        return Poly_ID

    if int(self_poly[6]) == 0:
        LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
        return Poly_ID
    else:
        if out_nbr_cnt > 2:
            LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
            return Poly_ID
        elif out_nbr_cnt == 1:
            LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
            LabelContLevel(inContour_Poly_ID,int(nbr_poly[0]),level_id)
            return int(nbr_poly[0])
        elif out_nbr_cnt == 0:
            LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
            return Poly_ID
        elif out_nbr_cnt == 2:
            LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
            out_nbr_nbr_ID = int(nbr_poly[7])
            if out_nbr_nbr_ID == 0:
                LabelContLevel(inContour_Poly_ID,Poly_ID,level_id)
                return Poly_ID
            else:
                return FindSameLevelUpperBound(inContour_Poly_ID,int(nbr_poly[0]),level_id)

def FindUpperLevel(inContour_Poly_ID,Poly_ID,level_id):
    self_poly = ReadPolyInfo(inContour_Poly_ID,Poly_ID)
    if int(self_poly[7]) > 0:
        nbr_poly = ReadPolyInfo(inContour_Poly_ID,int(self_poly[7]))
    else:
        return 0

    if int(nbr_poly[5]) < level_id and float(nbr_poly[1]) > float(self_poly[1]):
        LabelContLevel(inContour_Poly_ID,int(nbr_poly[0]),level_id)
        return int(nbr_poly[0])
    else:
        return 0



def GetContourID_Contour_Area(infc,areaThresholdValue):
    arcpy.AddField_management(infc, "Cont_ID", "LONG")
    arcpy.AddField_management(infc, "Cont_", "DOUBLE")
    arcpy.AddField_management(infc, "Lap_Area", "DOUBLE")

    feature_array = arcpy.Array()
    count = 0

    rows = arcpy.UpdateCursor(infc)
    for row in rows:

        # Create the geometry object
        feature = row.getValue("Shape")
        firstPnt = feature.firstPoint
        lastPnt = feature.lastPoint

        if  (firstPnt.X == lastPnt.X) & (firstPnt.Y == lastPnt.Y):
            count += 1
            #arcpy.AddMessage count
            feature_array.removeAll()
            for part in feature:
                for pnt in part:
                    if pnt:
                        feature_array.add(pnt)
                    else:
                    # If pnt is None, this represents an interior ring
                        arcpy.AddMessage("Interior Ring:")
            polygon = arcpy.Polygon(feature_array)
            area = polygon.area
            if(area < areaThresholdValue):
                area = 0
        else:
           area = 0
        row.setValue("Lap_Area",area)
        row.setValue("Cont_ID",row.getValue("FID"))
        row.setValue("Cont_",row.getValue("Contour"))
        rows.updateRow(row)
    outTempLayer = os.path.join(os.path.split(infc)[0],"tempLayer")
    arcpy.MakeFeatureLayer_management(infc,outTempLayer)
    selectExpression = arcpy.AddFieldDelimiters(outTempLayer,"Lap_Area") + "=0"
    arcpy.SelectLayerByAttribute_management(outTempLayer,"NEW_SELECTION",selectExpression)
    arcpy.DeleteFeatures_management(outTempLayer)

    arcpy.AddMessage("GetContourID_Contour_Area Over!")

def CalculateContour_Line_FromLine(outContours,gMinLength):
    arcpy.CheckOutExtension("Spatial")
    arcpy.AddField_management(outContours,"length","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management(outContours,"length","!shape.length@meters!","PYTHON_9.3","#")

    arcpy.AddField_management(outContours, "Distance", "DOUBLE")
    expression = "MySub( !Shape!)"
    codeblock = """def MySub(shape):
        fst = arcpy.PointGeometry(shape.firstPoint)
        lst = arcpy.PointGeometry(shape.lastPoint)
        return float(fst.distanceTo(lst))"""
    arcpy.CalculateField_management(outContours, "Distance", expression, "PYTHON_9.3", codeblock)


    tempLayer = os.path.join(os.path.split(outContours)[0],"temp_Line_Layer")
    arcpy.MakeFeatureLayer_management(outContours, tempLayer)
    expression = "length < " + str(gMinLength) + " or Distance > 0.01"
    # Execute SelectLayerByAttribute to determine which features to delete
    arcpy.SelectLayerByAttribute_management(tempLayer, "NEW_SELECTION",expression)

    # Execute GetCount and if some features have been selected, then
    #  execute DeleteFeatures to remove the selected features.
    if int(arcpy.GetCount_management(tempLayer).getOutput(0)) > 0:
        arcpy.DeleteFeatures_management(tempLayer)
    arcpy.Delete_management(tempLayer)
    arcpy.AddMessage("CalculateContour_Line Over!")
    #arcpy.AddMessage("CalculateContour_Line Over!")





def CalculateContour_Line(inRaster, outContours, contourInterval, baseContour,gMinLength):
    arcpy.CheckOutExtension("Spatial")
    arcpy.sa.Contour(inRaster, outContours, contourInterval, baseContour)
    #arcpy.sa.ContourWithBarriers(inRaster,outContours,"","POLYLINES","","",baseContour,contourInterval)
    arcpy.AddField_management(outContours,"length","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management(outContours,"length","!shape.length@meters!","PYTHON_9.3","#")

    arcpy.AddField_management(outContours, "Distance", "DOUBLE")
    expression = "MySub(!Shape!)"
    codeblock = """def MySub(shape):
        fst = arcpy.PointGeometry(shape.firstPoint)
        lst = arcpy.PointGeometry(shape.lastPoint)
        return float(fst.distanceTo(lst))"""
    arcpy.CalculateField_management(outContours, "Distance", expression, "PYTHON_9.3", codeblock)


    tempLayer = os.path.join(os.path.split(inRaster)[0],"temp_Line_Layer")
    arcpy.MakeFeatureLayer_management(outContours, tempLayer)
    expression = "length < " + str(gMinLength) + " or Distance > 0.01"
    # Execute SelectLayerByAttribute to determine which features to delete
    arcpy.SelectLayerByAttribute_management(tempLayer, "NEW_SELECTION",expression)

    # Execute GetCount and if some features have been selected, then
    #  execute DeleteFeatures to remove the selected features.
    if int(arcpy.GetCount_management(tempLayer).getOutput(0)) > 0:
        arcpy.DeleteFeatures_management(tempLayer)
    arcpy.Delete_management(tempLayer)
    arcpy.AddMessage("CalculateContour_Line Over!")
    #arcpy.AddMessage("CalculateContour_Line Over!")

def GetContourID_Contour_Area2(infc,outfc,areaThresholdValue):
    workspace = os.path.split(outfc)[0]
    outfcName = os.path.split(outfc)[1]

    feature_array = arcpy.Array()
    #count = 0
    origin_des = arcpy.Describe(infc)
    sr = origin_des.spatialReference
    arcpy.CreateFeatureclass_management(workspace,outfcName,"POLYGON","","","",sr)

    arcpy.AddField_management(infc, "Cont_ID", "LONG")
    arcpy.AddField_management(infc, "Cont_", "DOUBLE")
    arcpy.AddField_management(infc, "Lap_Area", "DOUBLE")

    arcpy.AddField_management(outfc, "Cont_ID", "LONG")
    arcpy.AddField_management(outfc, "Cont_", "DOUBLE")
    arcpy.AddField_management(outfc, "Lap_Area", "DOUBLE")
    # arcpy.AddField_management(outfc,"Lap_Area","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")

    cursor = arcpy.InsertCursor(outfc)
    rows = arcpy.UpdateCursor(infc)
    for row in rows:
        feature = row.getValue("Shape")
        Cont_ID = int(row.getValue("ID")) + 1
        feature_array.removeAll()
        for part in feature:
            for pnt in part:
                if pnt:
                    feature_array.add(pnt)

        f_polygon = arcpy.Polygon(feature_array,sr)

        area = f_polygon.area
        if (area >= float(areaThresholdValue)):
            newrow = cursor.newRow()
            newrow.SHAPE = f_polygon
            newrow.setValue("Cont_ID",Cont_ID)
            newrow.setValue("Cont_",row.getValue("Contour"))
            newrow.setValue("Lap_Area",area)
            cursor.insertRow(newrow)
        else:
            area = 0

        row.setValue("Cont_",row.getValue("Contour"))
        row.setValue("Lap_Area",area)
        row.setValue("Cont_ID",Cont_ID)
        rows.updateRow(row)
    del cursor

    outTempLayer = os.path.join(os.path.split(infc)[0],"tempLayer")
    arcpy.MakeFeatureLayer_management(infc,outTempLayer)
    selectExpression = "Lap_Area = 0 "
    arcpy.SelectLayerByAttribute_management(outTempLayer,"NEW_SELECTION",selectExpression)
    if int(arcpy.GetCount_management(outTempLayer).getOutput(0)) > 0:
        arcpy.DeleteFeatures_management(outTempLayer)
    arcpy.Delete_management(outTempLayer)
    arcpy.AddMessage("GetContourID_Contour_Area Over!")

def CalculateContour_Polygon(inContour_Line,outContour_Polygon,clusTol):
    arcpy.FeatureToPolygon_management(inContour_Line, outContour_Polygon,clusTol,"ATTRIBUTES", "")

    arcpy.AddMessage("CalculateContour_Polygon Over!")

def GetContourID4Polygon(inContourPolygon,inContourLapPolygon, outContourPolygonWithID,ToLevelID):
    if (os.path.splitext(inContourPolygon)[1].lower() == ".shp"):
        outContours_Polygon_WithID_Temp = os.path.join(os.path.split(outContourPolygonWithID)[0] ,"IDpolygonTemp.shp")
    else:
        outContours_Polygon_WithID_Temp = os.path.join(os.path.split(outContourPolygonWithID)[0] ,"IDpolygonTemp")
    fieldmappings = arcpy.FieldMappings()
    fieldmappings.addTable(inContourLapPolygon)
    fieldmappings.addTable(inContourPolygon)
    arcpy.SpatialJoin_analysis(inContourPolygon,inContourLapPolygon,outContours_Polygon_WithID_Temp,"JOIN_ONE_TO_ONE","KEEP_ALL",fieldmappings,"WITHIN_CLEMENTINI")
    arcpy.Dissolve_management(outContours_Polygon_WithID_Temp,outContourPolygonWithID,["Cont_ID"],[["Cont_","FIRST"],["Lap_Area","SUM"]])
    arcpy.Delete_management(outContours_Polygon_WithID_Temp)

    arcpy.AddField_management(outContourPolygonWithID, "Line_Cou", "LONG")
    arcpy.AddField_management(outContourPolygonWithID,"Level_ID","SHORT")
    arcpy.AddField_management(outContourPolygonWithID,"Level_BD","SHORT")
    arcpy.AddField_management(outContourPolygonWithID,"LabelFlag","SHORT")
    arcpy.AddField_management(outContourPolygonWithID,"Out_Nbr_Ct","SHORT")
    arcpy.AddField_management(outContourPolygonWithID,"Out_Nbr_ID","SHORT")
    arcpy.AddField_management(outContourPolygonWithID,"Contain_ID","TEXT")

    expression = "MySub( !Shape!)"
    codeblock = """def MySub(feat):
        if feat.isMultipart:
            return 2
        else:
            return 1"""
    arcpy.CalculateField_management(outContourPolygonWithID, "Line_Cou", expression, "PYTHON_9.3", codeblock)
    arcpy.AddMessage("GetContourID4Polygon Over!")

def AlterField_management(inFeatureClass, outFeatureClass, new_name_by_old_name):
    """ Renames specified fields in input feature class/table
    :table:                 input table (fc, table, layer, etc)
    :out_table:             output table (fc, table, layer, etc)
    :new_name_by_old_name:  {'old_field_name':'new_field_name',...}
    ->  out_table
    """
    existing_field_names = [field.name for field in arcpy.ListFields(inFeatureClass)]

    field_mappings = arcpy.FieldMappings()
    field_mappings.addTable(inFeatureClass)

    for old_field_name, new_field_name in new_name_by_old_name.iteritems():
        if old_field_name not in existing_field_names:
            message = "Field: {0} not in {1}".format(old_field_name, inFeatureClass)
            raise Exception(message)

        mapping_index = field_mappings.findFieldMapIndex(old_field_name)
        field_map = field_mappings.fieldMappings[mapping_index]
        output_field = field_map.outputField
        output_field.name = new_field_name
        output_field.aliasName = new_field_name
        field_map.outputField = output_field
        field_mappings.replaceFieldMap(mapping_index, field_map)

    # use merge with single input just to use new field_mappings
    arcpy.Merge_management(inFeatureClass, outFeatureClass, field_mappings)
    return outFeatureClass

def GetPolygonNeighbor(inContour_Polygon_ID, outContours_Polygon_Neigh,outContours_Polygon_Summarize):
    workSpace = os.path.split(outContours_Polygon_Neigh)[0]
    if (os.path.splitext(outContours_Polygon_Neigh)[1].lower() == ".dbf"):
        outContours_Polygon_Neigh_Temp = os.path.join(workSpace,"tempNei.dbf")
        arcpy.PolygonNeighbors_analysis(inContour_Polygon_ID,outContours_Polygon_Neigh_Temp,"FID;Cont_ID;FIRST_Cont;SUM_Lap_Ar", "NO_AREA_OVERLAP", "BOTH_SIDES")
        arcpy.DeleteField_management(outContours_Polygon_Neigh_Temp,["src_FID","nbr_FID","LENGTH","NODE_COUNT"])
        new_name_by_old_name = { 'src_FIRST_':'src_Cont','nbr_FIRST_':'nbr_Cont','src_SUM_La':'src_Area','nbr_SUM_La':'nbr_Area' }
    else:
        outContours_Polygon_Neigh_Temp = os.path.join(workSpace,"tempNei")
        arcpy.PolygonNeighbors_analysis(inContour_Polygon_ID,outContours_Polygon_Neigh_Temp,"Cont_ID;FIRST_Cont_;SUM_Lap_Area", "NO_AREA_OVERLAP", "BOTH_SIDES")
        arcpy.DeleteField_management(outContours_Polygon_Neigh_Temp,["LENGTH","NODE_COUNT"])
        new_name_by_old_name = { 'src_FIRST_Cont_':'src_Cont','nbr_FIRST_Cont_':'nbr_Cont','src_SUM_Lap_Area':'src_Area','nbr_SUM_Lap_Area':'nbr_Area' }

    AlterField_management(outContours_Polygon_Neigh_Temp,outContours_Polygon_Neigh,new_name_by_old_name)
    arcpy.Delete_management(outContours_Polygon_Neigh_Temp)

    workSpace2 = os.path.split(outContours_Polygon_Summarize)[0]
    if (os.path.splitext(outContours_Polygon_Summarize)[1].lower() == ".dbf"):
        outContours_Polygon_Summarize_Temp = os.path.join(workSpace2,"tempSum.dbf")
        arcpy.Statistics_analysis(outContours_Polygon_Neigh,outContours_Polygon_Summarize_Temp,[["src_Cont_I", "COUNT"]],"src_Cont_I")
        arcpy.DeleteField_management(outContours_Polygon_Summarize_Temp,"COUNT_src_")
    else:
        outContours_Polygon_Summarize_Temp = os.path.join(workSpace2,"tempSum")
        arcpy.Statistics_analysis(outContours_Polygon_Neigh,outContours_Polygon_Summarize_Temp,[["src_Cont_ID", "COUNT"]],"src_Cont_ID")
        arcpy.DeleteField_management(outContours_Polygon_Summarize_Temp,"COUNT_src_Cont_ID")

    AlterField_management(outContours_Polygon_Summarize_Temp,outContours_Polygon_Summarize,{'FREQUENCY':'NBR_CNT'})
    arcpy.Delete_management(outContours_Polygon_Summarize_Temp)
    arcpy.AddMessage("GetPolygonNeighbor Over!")

def GetStatisticalDataTable(inDEMRaster,inContoutPolygonID,outStatisticalTable):
    arcpy.CheckOutExtension("Spatial")

    raster = arcpy.Raster(inDEMRaster)
    cellSizeHeight = raster.meanCellHeight

    if (os.path.splitext(inContoutPolygonID)[1].lower() == ".shp"):
        outZoneRaster = os.path.join(os.path.split(inContoutPolygonID)[0],"zoneRaster2.tif")
    else:
        database_space = os.path.split(os.path.split(inContoutPolygonID)[0])
        outZoneRaster = os.path.join(os.path.split(inContoutPolygonID)[0],"zoneRaster2")
    arcpy.CheckOutExtension("Spatial")
    #arcpy.sa.ZonalStatisticsAsTable(inContour_Poly_ID,"Cont_ID",inDEMRaster,outStatisticalTable,"DATA","ALL")
    arcpy.PolygonToRaster_conversion(inContoutPolygonID,"Cont_ID",outZoneRaster,"CELL_CENTER","",cellSizeHeight)

    arcpy.sa.ZonalStatisticsAsTable(outZoneRaster,"Value",inDEMRaster,outStatisticalTable,"DATA","ALL")
    arcpy.AddField_management(outStatisticalTable, "Cont_ID", "LONG")
    arcpy.CalculateField_management(outStatisticalTable,"Cont_ID","!Value!","PYTHON_9.3")
    arcpy.Delete_management(outZoneRaster)

    arcpy.AddMessage("GetStatisticalDataTable Over!")

def numpyToArray(AllInfo,inContour_Poly_ID):
    if (os.path.splitext(inContour_Poly_ID)[1].lower() == ".shp"):
        outTempTable = os.path.join(os.path.split(inContour_Poly_ID)[0],"tempTable.dbf")
    else:
        outTempTable = os.path.join(os.path.split(inContour_Poly_ID)[0],"tempTable")
    arcpy.da.NumPyArrayToTable(AllInfo, outTempTable)
    arcpy.JoinField_management(inContour_Poly_ID,"Cont_ID",outTempTable,"Cont_ID")
    arcpy.CalculateField_management(inContour_Poly_ID,"Level_ID","!Level_ID_1!","PYTHON_9.3")
    arcpy.CalculateField_management(inContour_Poly_ID,"Level_BD","!Level_BD_1!","PYTHON_9.3")
    #arcpy.RemoveJoin_management(inContour_Poly_ID, outTempTable)
    arcpy.Delete_management(outTempTable)


def LapJoin(inContour_Polygon_ID,inContour_Polygon_Lap_temp):
    # arcpy.Dissolve_management(inContour_Polygon_Lap,inContour_Polygon_Lap_temp,["Cont_ID"],[["Cont_","FIRST"],["Lap_Area","SUM"]])
    # arcpy.Delete_management(inContour_Polygon_Lap)
    arcpy.AddIndex_management (inContour_Polygon_Lap_temp, "Cont_ID", "Cont_Index")
    arcpy.JoinField_management(inContour_Polygon_Lap_temp,"Cont_ID",inContour_Polygon_ID,"Cont_ID",["Level_BD"])


def SingleShapefile(inContour_Polygon_Lap_temp,OutPath,Max_level,inDEMRaster,isContain):
    if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
        arcpy.CreateFolder_management(OutPath, "Single")
    else:
        sr = arcpy.Describe(inContour_Polygon_Lap_temp).spatialReference
        arcpy.CreateFeatureDataset_management(OutPath, "Single",sr)

    Single_OutPath = os.path.join(OutPath,"Single")
    for level_id in range(1,int(Max_level)):
        selection = "Level_BD = " + str(level_id)
        if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
            shapefileName = "Single_L" + str(level_id) + ".shp"
            outputTemp = os.path.join(Single_OutPath,"temp_" + str(level_id) + ".shp")
        else:
            shapefileName = "Single_L" + str(level_id)
            outputTemp = os.path.join(Single_OutPath,"temp_" + str(level_id))
        outputShapefile = os.path.join(Single_OutPath,shapefileName)

        arcpy.MakeFeatureLayer_management (inContour_Polygon_Lap_temp,outputTemp )
        arcpy.SelectLayerByAttribute_management (outputTemp, "NEW_SELECTION", selection)
        arcpy.CopyFeatures_management(outputTemp,outputShapefile)
        arcpy.Delete_management(outputTemp)

    for level_id in range(1,int(Max_level)):
        if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
            shapefileName = "Single_L" + str(level_id) + ".shp"
        else:
            shapefileName = "Single_L" + str(level_id)
        outputShapefile = os.path.join(Single_OutPath,shapefileName)
        temp_level_id = level_id
        while(temp_level_id < int(Max_level)-1):
            temp_level_id += 1
            if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
                shapefileNameTemp = "Single_L" + str(temp_level_id) + ".shp"
                out_feature_class  = os.path.join(Single_OutPath,"SingleSJTemp.shp")
            else:
                shapefileNameTemp = "Single_L" + str(temp_level_id)
                out_feature_class  = os.path.join(Single_OutPath,"SingleSJTemp")

            outputShapefileTemp = os.path.join(Single_OutPath,shapefileNameTemp)


            arcpy.SpatialJoin_analysis(outputShapefileTemp,outputShapefile,out_feature_class, "JOIN_ONE_TO_MANY","KEEP_COMMON",match_option="WITHIN")
            result = arcpy.GetCount_management(out_feature_class)
            overlayCount = int(result.getOutput(0))
            if overlayCount>0:
                if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
                    outEraseTemp = os.path.join(Single_OutPath,"EraseTemp.shp")
                else:
                    outEraseTemp = os.path.join(Single_OutPath,"EraseTemp")
                arcpy.Erase_analysis(outputShapefile,out_feature_class,outEraseTemp)
                arcpy.Delete_management(outputShapefile)
                arcpy.Rename_management(outEraseTemp,outputShapefile)
            arcpy.Delete_management(out_feature_class)


        arcpy.CalculateField_management(outputShapefile,"Lap_Area","!shape.area@squaremeters!","PYTHON_9.3","#")
        CalculateVolume(outputShapefile,inDEMRaster)

        if isContain:
            Contain_OutPath = os.path.join(OutPath,"Contains")
            if level_id == 1:
                if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
                    arcpy.CreateFolder_management(OutPath, "Contains")
                    contain_ShapefileName = "Contain_L" + str(level_id) + ".shp"
                else:
                    sr = arcpy.Describe(inContour_Polygon_Lap_temp).spatialReference
                    arcpy.CreateFeatureDataset_management(OutPath, "Contains",sr)
                    contain_ShapefileName = "Contain_L" + str(level_id)

                outputContainShapefile = os.path.join(Contain_OutPath,contain_ShapefileName)
                arcpy.CopyFeatures_management(outputShapefile,outputContainShapefile)
            else:
                if (os.path.splitext(inContour_Polygon_Lap_temp)[1].lower() == ".shp"):
                    contain_ShapefileName = "Contain_L" + str(level_id) + ".shp"
                    preShape = os.path.join(Contain_OutPath,"Contain_L" + str(level_id-1) + ".shp")
                else:
                    contain_ShapefileName = "Contain_L" + str(level_id)
                    preShape = os.path.join(Contain_OutPath,"Contain_L" + str(level_id-1))
                outputContainShapefile = os.path.join(Contain_OutPath,contain_ShapefileName)

                arcpy.Update_analysis(preShape,outputShapefile,outputContainShapefile)
                arcpy.CalculateField_management(outputContainShapefile,"Lap_Area","!shape.area@squaremeters!","PYTHON_9.3","#")
                CalculateVolume(outputContainShapefile,inDEMRaster)

def IndividualShapefile(inContour_Polygon_Lap_temp,OutPath,Max_level,isContain):
    flag = True
    for level_id in range(1,int(Max_level)):
        if(isContain):
            if flag:
                arcpy.CreateFolder_management(OutPath, "Contains")
                OutPath = os.path.join(OutPath,"Contains")
                flag = False
            selection = "Level_BD <= " + str(level_id) + " and Level_BD > 0"
            shapefileName = "Contains_L" + str(level_id) + ".shp"
        else:
            if flag:
                arcpy.CreateFolder_management(OutPath, "Single")
                OutPath = os.path.join(OutPath,"Single")
                flag = False
            selection = "Level_BD = " + str(level_id)
            shapefileName = "Single_L" + str(level_id) + ".shp"

        outputShapefile = os.path.join(OutPath,shapefileName)
        outputTemp = os.path.join(OutPath,"temp_" + str(level_id) + ".shp")
        arcpy.MakeFeatureLayer_management (inContour_Polygon_Lap_temp,outputTemp )
        arcpy.SelectLayerByAttribute_management (outputTemp, "NEW_SELECTION", selection)
        arcpy.CopyFeatures_management(outputTemp,outputShapefile)
        arcpy.Delete_management(outputTemp)


def CreateTIN(inLine,outTin):
    arcpy.CheckOutExtension("3D")
    arcpy.CreateTin_3d(outTin,in_features=inLine)




#Parameters
#Input Data
DEMRasterPath = arcpy.GetParameterAsText(0)
contourInterval = float(arcpy.GetParameterAsText(1))
baseContour = float(arcpy.GetParameterAsText(2))
gMinArea = float(arcpy.GetParameterAsText(3))
gMinDepth = float(arcpy.GetParameterAsText(4))
OutputPath = arcpy.GetParameterAsText(5)
gMinLength = 40
#min_size = 100
#buff_dis = 5
#OutputDB = r"G:\Chenzuoqi\ContourTree\DataTest\ContourTemp"
workspace =os.path.split(OutputPath)[0]
isContain = True
doPeakProcess = True

#Temp Data and Output Data
inDEMRaster = DEMRasterPath

if (os.path.splitext(OutputPath)[1].lower() == ".shp"):
    outContours_Line = os.path.join(workspace,"outContours_Line.shp")
    outContours_Polygon = os.path.join(workspace,"outContours_Polygon.shp")
    clusTol = "0"
    outContours_Polygon_WithID = OutputPath
    outContours_Polygon_Neigh = os.path.join(workspace,"outPolygon_Nei.dbf")
    outContours_Polygon_Summarize = os.path.join(workspace,"outPolygon_Sum.dbf")
    outContours_Polygon_Overlap = os.path.join(workspace,"outPolygon_Lap.shp")
    outContours_Statistical_Table = os.path.join(workspace,"outStatistic.dbf")
else:
    outContours_Line = os.path.join(workspace,"outContours_Line")
    outContours_Polygon = os.path.join(workspace,"outContours_Polygon")
    clusTol = "0"
    outContours_Polygon_WithID = OutputPath
    outContours_Polygon_Neigh = os.path.join(workspace,"outPolygon_Nei")
    outContours_Polygon_Summarize = os.path.join(workspace,"outPolygon_Sum")
    outContours_Polygon_Overlap = os.path.join(workspace,"outPolygon_Lap")
    outContours_Statistical_Table = os.path.join(workspace,"outStatistic")
ToLevelID = []
HasLevelID = []


beginDT = datetime.datetime.now()
beginDTAll = beginDT
print ("Begin, Now is %s}" %(datetime.datetime.strftime(beginDT,"%Y-%m-%d %H:%M:%S")))

#inDEMRaster = PreProcess(DEMRasterPath,min_size,buff_dis)

CalculateContour_Line(inDEMRaster,outContours_Line,contourInterval,baseContour,gMinLength)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

GetContourID_Contour_Area2(outContours_Line,outContours_Polygon_Overlap,gMinArea)
# #GetContourID_Contour_Area(outContours_Line,areaThresholdValue)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

CalculateContour_Polygon(outContours_Line,outContours_Polygon,clusTol)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

GetContourID4Polygon(outContours_Polygon,outContours_Polygon_Overlap, outContours_Polygon_WithID,ToLevelID)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

GetPolygonNeighbor(outContours_Polygon_WithID, outContours_Polygon_Neigh,outContours_Polygon_Summarize)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

GetStatisticalDataTable(inDEMRaster,outContours_Polygon_WithID,outContours_Statistical_Table)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

JoinAllTables(outContours_Polygon_WithID,outContours_Statistical_Table,outContours_Polygon_Summarize,outContours_Line)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT


AllInfo = ReadAllInfo(outContours_Polygon_WithID)

UpdateAllNbrInfo(AllInfo,outContours_Polygon_Neigh)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

IdentifyLevel0(AllInfo)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

IdentifyLevel1(AllInfo,gMinDepth,gMinArea)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

IdentifyLevel2(AllInfo)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT

Max_level = IdentifyAllLevel(AllInfo)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))
beginDT = endDT


if doPeakProcess:
    IdentifyPeakLevel(AllInfo)
    endDT = datetime.datetime.now()
    deltaDT = (endDT - beginDT).seconds
    print ("Using %s s" %(deltaDT))
    beginDT = endDT

PostProcessing(AllInfo,gMinDepth,gMinArea)
endDT = datetime.datetime.now()
deltaDT = (endDT - beginDT).seconds
print ("Using %s s" %(deltaDT))

numpyToArray(AllInfo,outContours_Polygon_WithID)
LapJoin(outContours_Polygon_WithID,outContours_Polygon_Overlap)
SingleShapefile(outContours_Polygon_Overlap,workspace,Max_level,inDEMRaster,isContain)
deltaDT = (endDT - beginDTAll).seconds
print ("Total Using %s s" %(deltaDT))



if not (os.path.splitext(outContours_Polygon)[1].lower() == ".shp"):
    arcpy.Compact_management(workspace)
arcpy.Delete_management(outContours_Line)
arcpy.Delete_management(outContours_Polygon)
arcpy.Delete_management(outContours_Polygon_Neigh)
arcpy.Delete_management(outContours_Polygon_Summarize)
arcpy.Delete_management(outContours_Polygon_Overlap)
arcpy.Delete_management(outContours_Statistical_Table)
print ("All Over!")
