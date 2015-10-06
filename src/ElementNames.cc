#include "../include/ElementNames.hh"

string* ElementNames::elementName=NULL;

using namespace std;

ElementNames::ElementNames()
{
    //elementName=NULL;
}

ElementNames::~ElementNames()
{
    //dtor
}

void ElementNames::ClearStore()
{
    if(elementName != NULL)
        delete [] elementName;
}

void ElementNames::SetElementNames()
{
    elementName = new string[119];

    elementName[0] = "Error";
    elementName[1] = "Hydrogen";
    elementName[2] = "Helium";
    elementName[3] = "Lithium";
    elementName[4] = "Beryllium";
    elementName[5] = "Boron";
    elementName[6] = "Carbon";
    elementName[7] = "Nitrogen";
    elementName[8] = "Oxygen";
    elementName[9] = "Fluorine";
    elementName[10] = "Neon";
    elementName[11] = "Sodium";
    elementName[12] = "Magnesium";
    elementName[13] = "Aluminum";
    elementName[14] = "Silicon";
    elementName[15] = "Phosphorous";
    elementName[16] = "Sulfur";
    elementName[17] = "Chlorine";
    elementName[18] = "Argon";
    elementName[19] = "Potassium";
    elementName[20] = "Calcium";
    elementName[21] = "Scandium";
    elementName[22] = "Titanium";
    elementName[23] = "Vanadium";
    elementName[24] = "Chromium";
    elementName[25] = "Manganese";
    elementName[26] = "Iron";
    elementName[27] = "Cobalt";
    elementName[28] = "Nickel";
    elementName[29] = "Copper";
    elementName[30] = "Zinc";
    elementName[31] = "Gallium";
    elementName[32] = "Germanium";
    elementName[33] = "Arsenic";
    elementName[34] = "Selenium";
    elementName[35] = "Bromine";
    elementName[36] = "Krypton";
    elementName[37] = "Rubidium";
    elementName[38] = "Strontium";
    elementName[39] = "Yttrium";
    elementName[40] = "Zirconium";
    elementName[41] = "Niobium";
    elementName[42] = "Molybdenum";
    elementName[43] = "Technetium";
    elementName[44] = "Ruthenium";
    elementName[45] = "Rhodium";
    elementName[46] = "Palladium";
    elementName[47] = "Silver";
    elementName[48] = "Cadmium";
    elementName[49] = "Indium";
    elementName[50] = "Tin";
    elementName[51] = "Antimony";
    elementName[52] = "Tellurium";
    elementName[53] = "Iodine";
    elementName[54] = "Xenon";
    elementName[55] = "Cesium";
    elementName[56] = "Barium";
    elementName[57] = "Lanthanum";
    elementName[58] = "Cerium";
    elementName[59] = "Praseodymium";
    elementName[60] = "Neodymium";
    elementName[61] = "Promethium";
    elementName[62] = "Samarium";
    elementName[63] = "Europium";
    elementName[64] = "Gadolinium";
    elementName[65] = "Terbium";
    elementName[66] = "Dysprosium";
    elementName[67] = "Holmium";
    elementName[68] = "Erbium";
    elementName[69] = "Thulium";
    elementName[70] = "Ytterbium";
    elementName[71] = "Lutetium";
    elementName[72] = "Hafnium";
    elementName[73] = "Tantalum";
    elementName[74] = "Tungsten";
    elementName[75] = "Rhenium";
    elementName[76] = "Osmium";
    elementName[77] = "Iridium";
    elementName[78] = "Platinum";
    elementName[79] = "Gold";
    elementName[80] = "Mercury";
    elementName[81] = "Thallium";
    elementName[82] = "Lead";
    elementName[83] = "Bismuth";
    elementName[84] = "Polonium";
    elementName[85] = "Astatine";
    elementName[86] = "Radon";
    elementName[87] = "Francium";
    elementName[88] = "Radium";
    elementName[89] = "Actinium";
    elementName[90] = "Thorium";
    elementName[91] = "Protactinium";
    elementName[92] = "Uranium";
    elementName[93] = "Neptunium";
    elementName[94] = "Plutonium";
    elementName[95] = "Americium";
    elementName[96] = "Curium";
    elementName[97] = "Berkelium";
    elementName[98] = "Californium";
    elementName[99] = "Einsteinium";
    elementName[100] = "Fermium";
    elementName[101] = "Mendelevium";
    elementName[102] = "Nobelium";
    elementName[103] = "Lawrencium";
    elementName[104] = "Rutherfordium";
    elementName[105] = "Dubnium";
    elementName[106] = "Seaborgium";
    elementName[107] = "Bohrium";
    elementName[108] = "Hassium";
    elementName[109] = "Meitnerium";
    elementName[110] = "Darmstadtium";
    elementName[111] = "Roentgenium";
    elementName[112] = "Copernicium";
    elementName[113] = "Ununtrium";
    elementName[114] = "Flerovium";
    elementName[115] = "Ununpentium";
    elementName[116] = "Livermorium";
    elementName[117] = "Ununseptium";
    elementName[118] = "Ununoctium";

}

bool ElementNames::CheckName(string name, int Z)
{

    if(elementName != NULL)
    {
        if(name.substr((name.length()-2),2)==".z")
        {
            name.erase((name.length()-2),2);
        }

        if(Z==0)
            return false;

        if(name == elementName[Z])
            return true;
        else
        {
            if(StringPercentMatch(elementName[Z],name)>=0.8)
                return true;
        }
    }
    else
    {
        cout << "### Error SetElementNames() has not been yet ###" << endl;
    }


    return false;

}

bool ElementNames::CheckName(string name)
{
    if(elementName != NULL)
    {
        if(name.substr((name.length()-2),2)==".z")
        {
            name.erase((name.length()-2),2);
        }

        for(int i=0; i<119; i++)
        {
            if(name == elementName[i])
                return true;
            else
            {
                name[0] = name[0]+('A'-'a');
                if(name == elementName[i])
                    return true;
                else if(StringPercentMatch(elementName[i],name)>=0.8)
                    return true;
            }
        }
    }
    else
    {
        cout << "### Error SetElementNames() has not been yet ###" << endl;
    }

    return false;

}

double ElementNames::StringPercentMatch(string base, string secondary)
{
    int size1, size2, match, sizeB, pos1, pos2, count1=-1, count2=-1, tally2=0;
    string nameMatch=secondary, steps="";
    stringstream temp;
    double percent, sum, tally=0;

    temp.str(base);
    size1=(ExtractString(temp,'\n',4)).size()-1;
    if(!temp)
    {
        temp.clear();
    }
    temp.str(secondary);
    size2=(ExtractString(temp,'\n',4)).size()-1;

    sum=max(base.size()-0.5*size1, secondary.size()-0.5*size2);
    sizeB=max(base.size(), secondary.size());
    match=min(base.size(), secondary.size());

    for(int l=0; l<sizeB; l++)
    {
        //don't adjust match here since it is already adjusted above
        if(l==int(nameMatch.size()))
        {
            nameMatch.push_back(base[l]);
            steps.push_back('+');
        }
        else if(nameMatch[l]!=base[l])
        {
            pos1=nameMatch.find_first_of(base[l],l);
            if(pos1!=int(string::npos))
            {
                count1=0;
                while(((l+count1)<sizeB)&&(pos1<int(nameMatch.size()))&&(nameMatch[pos1]==base[l+count1]))
                {
                    pos1++;
                    count1++;
                }
                pos1-=count1;
            }
            pos2=base.find_first_of(nameMatch[l],l);
            if(pos2!=int(string::npos))
            {
                count2=0;
                while((pos2<sizeB)&&((l+count2)<int(nameMatch.size()))&&(nameMatch[l+count2]==base[pos2]))
                {
                    pos2++;
                    count2++;
                }
                pos2-=count2;
            }

            if((pos1!=int(string::npos))&&(pos2!=int(string::npos)))
            {
                if(((count1-(pos1-l))>=(count2-(pos2-l)))&&((count1-(pos1-l))>=0))
                {
                    nameMatch.erase(l,pos1-l);
                    steps.append(pos1-l, '-');
                    tally2+=pos1-l;
                    match-=pos1-l;
                    l--;
                }
                else if((count2-(pos2-l))>=0)
                {
                    nameMatch.insert(l,base, l, pos2-l);
                    tally2+=pos2-l;
                    steps.append(pos2-l, '+');
                    match-=pos2-l;
                    l+=pos2-l-1;
                }
                else
                {
                    nameMatch[l]=base[l];
                    tally2++;
                    steps.push_back('s');
                    match--;
                }
            }
            else if(pos2!=int(string::npos))
            {
                if((count2-(pos2-l))>=0)
                {
                    nameMatch.insert(l,base, l, pos2-l);
                    tally2+=pos2-l;
                    steps.append(pos2-l, '+');
                    match-=pos2-l;
                    l+=pos2-l-1;
                }
                else
                {
                    nameMatch[l]=base[l];
                    tally2++;
                    steps.push_back('s');
                    match--;
                }
            }
            else if(pos1!=int(string::npos))
            {
                if((count1-(pos1-l))>=0)
                {
                    nameMatch.erase(l,pos1-l);
                    tally2+=pos1-l;
                    steps.append(pos1-l, '-');
                    match-=pos1-l;
                    l--;
                }
                else
                {
                    nameMatch[l]=base[l];
                    tally2++;
                    steps.push_back('s');
                    match--;
                }
            }
            else
            {
                nameMatch[l]=base[l];
                tally2++;
                steps.push_back('s');
                match--;
            }
        }
        else /*if(int(steps.size())!=int(nameMatch.size()))*/
        {
            steps.push_back('=');
            tally++;
            if(!(((nameMatch[l]>='A')&&(nameMatch[l]<='Z'))||((nameMatch[l]>='a')&&(nameMatch[l]<='z'))||((nameMatch[l]>='0')&&(nameMatch[l]<='9'))||(nameMatch[l]=='.')||(nameMatch[l]=='-')))
            {
                tally-=0.5;
            }
        }
    }

    if(int(nameMatch.size())>sizeB)
    {
        nameMatch.erase(sizeB,int(nameMatch.size())-sizeB);
        steps.append(int(nameMatch.size())-sizeB, '-');
    }

    percent = tally/sum;

    //cout << "\nbase: " << base << " secondary: " << secondary << " match: " << match << " tally " << tally << " tally2: " << tally2 << " base size: " << sizeB << " secondary size: " << secondary.size() << " steps: " << steps << endl;

    return percent;
}

string ElementNames::ExtractString(std::stringstream &stream, char delim, int outType)
{
    string value="";
    bool charOut=false, numOut=false, symOut=false;
    char letter;
    //bool first=true;

    if(outType==0)
    {

    }
    else if(outType==1)
    {
        charOut=true;
    }
    else if(outType==2)
    {
        numOut=true;
    }
    else if(outType==3)
    {
        charOut=true;
        numOut=true;
    }
    else if(outType==4)
    {
        symOut=true;
    }
    else if(outType==5)
    {
        charOut=true;
        symOut=true;
    }
    else if(outType==6)
    {
        numOut=true;
        symOut=true;
    }
    else
    {
        charOut=true;
        numOut=true;
        symOut=true;
    }

    while(stream&&(stream.peek()!=delim))
    {
        letter = stream.get();
        if(((letter>='A')&&(letter<='Z'))||((letter>='a')&&(letter<='z')))
        {
            if(charOut)
            {
                value+=letter;
                //first=true;
            }
            /*else if(first)
            {
                value+=' ';
                first=false;
            }*/
        }
        else if(((letter>='0')&&(letter<='9'))||(letter=='.')||(letter=='-'))
        {
            if(numOut)
            {
                value+=letter;
                //first=true;
            }
            /*else if(first)
            {
                value+=' ';
                first=false;
            }*/
        }
        else
        {
            if(symOut)
            {
                value+=letter;
                //first=true;
            }
            /*else if(first)
            {
                value+=' ';
                first=false;
            }*/
        }
    }
    return value;
}
