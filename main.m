function ans = main()
ans = 0;
continue_p = true;
    while continue_p == true
        filename = input("Input the name of file that you want to run: EchantillonageSn, ConvertCSVTOMAP, Question11, Saison \n\n"); 
        switch filename
            case "EchantillonageSn" 
               EchantillonageSn();
            case "Question11"
                Question11();
            case "Saison"
                Saison()
            case "ConvertCSVTOMAP"
                ConvertCSVTOMAP()
        end
        prompt = "Do you want more? Y/N [Y]: ";
        txt = input(prompt,"s");
        if isempty(txt)
            txt = 'Y';
        end
        continue_p = txt == 'Y';
    end
ans = 0;
exit(0);
end